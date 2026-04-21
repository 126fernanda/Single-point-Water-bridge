import json
import logging
import csv
import sys
import MDAnalysis as mda
import numpy as np

from .core import build_graph, compute_edge_probabilities, traverse_network

logger = logging.getLogger(__name__)

def cluster_pathways(data_file, threshold=3.5, output_file="clustered_pathways.json"):
    """
    Reads JSONL trajectory data and performs temporal clustering to identify collective pathways.
    """
    logger.info("Starting temporal clustering of pathways...")
    import scipy.spatial.distance as ssd
    from scipy.cluster.hierarchy import linkage, fcluster

    paths = []
    total_frames = 0
    frame_set = set()

    # Read paths into memory for clustering
    with open(data_file, 'r') as f:
        for line in f:
            obj = json.loads(line)
            if obj.get('type') == 'metadata':
                total_frames = obj.get('n_frames_analyzed', 0)
            elif obj.get('type') == 'frame':
                frame_idx = obj['frame_idx']
                frame_set.add(frame_idx)
                for p in obj['paths']:
                    paths.append({
                        'frame': frame_idx,
                        'nodes': p['nodes'],
                        'coords': np.array(p['coords']),
                        'probability': p['probability'],
                        'length': p['length']
                    })

    n_paths = len(paths)
    if total_frames == 0:
        total_frames = len(frame_set)

    if n_paths == 0:
        logger.info("No paths found to cluster.")
        return

    if n_paths == 1:
        logger.info("Only 1 path found. Skipping clustering.")
        return

    logger.info(f"Computing pairwise directed Hausdorff distance matrix for {n_paths} pathways...")

    # Compute condensed distance matrix
    # Note: directed_hausdorff returns (dist, index1, index2). We just need the distance.
    # Because lengths differ, Hausdorff is robust. We use the max of the two directed distances to make it symmetric.
    dist_matrix = np.zeros(n_paths * (n_paths - 1) // 2)
    idx = 0
    for i in range(n_paths):
        for j in range(i + 1, n_paths):
            u_coords = paths[i]['coords']
            v_coords = paths[j]['coords']
            d1 = ssd.directed_hausdorff(u_coords, v_coords)[0]
            d2 = ssd.directed_hausdorff(v_coords, u_coords)[0]
            dist_matrix[idx] = max(d1, d2)
            idx += 1

    logger.info("Performing hierarchical average-link clustering...")
    Z = linkage(dist_matrix, method='average')
    labels = fcluster(Z, t=threshold, criterion='distance')

    unique_labels = set(labels)
    n_clusters = len(unique_labels)
    logger.info(f"Identified {n_clusters} unique collective pathways.")

    clusters_data = []

    for label in unique_labels:
        cluster_indices = np.where(labels == label)[0]
        cluster_paths = [paths[i] for i in cluster_indices]

        # Calculate cluster occupancy
        cluster_frames = set(p['frame'] for p in cluster_paths)
        occupancy = len(cluster_frames) / total_frames if total_frames > 0 else 0.0

        # Calculate average probability
        avg_prob = np.mean([p['probability'] for p in cluster_paths])

        # Find the medoid (path with minimum average distance to all other paths in cluster)
        # For simplicity and speed if cluster is large, we can just pick the first one, or do a small inner loop
        if len(cluster_indices) <= 2:
            medoid_coords = cluster_paths[0]['coords'].tolist()
        else:
            min_dist_sum = float('inf')
            medoid_idx = 0
            for i, p1 in enumerate(cluster_paths):
                dist_sum = 0
                for p2 in cluster_paths:
                    d1 = ssd.directed_hausdorff(p1['coords'], p2['coords'])[0]
                    d2 = ssd.directed_hausdorff(p2['coords'], p1['coords'])[0]
                    dist_sum += max(d1, d2)
                if dist_sum < min_dist_sum:
                    min_dist_sum = dist_sum
                    medoid_idx = i
            medoid_coords = cluster_paths[medoid_idx]['coords'].tolist()

        clusters_data.append({
            "cluster_id": int(label),
            "size": len(cluster_indices),
            "occupancy": float(occupancy),
            "avg_probability": float(avg_prob),
            "medoid_coords": medoid_coords
        })

    # Sort by occupancy descending
    clusters_data.sort(key=lambda x: x['occupancy'], reverse=True)

    with open(output_file, 'w') as f:
        json.dump(clusters_data, f, indent=2)

    logger.info(f"Clustering complete. Results saved to {output_file}")


def sanitize_csv_field(field_value):
    """
    Sanitizes a field value for CSV export to prevent formula injection.
    If the value starts with =, +, -, or @, it prepends a single quote.
    """
    field_str = str(field_value)
    if field_str.startswith(('=', '+', '-', '@')):
        return f"'{field_str}"
    return field_str

def run_analysis(topo_file, traj_file, root_sel, water_sel="resname SOL or resname WAT or resname HOH",
                 stride=1, max_depth=10, min_depth=1, prob_threshold=1e-3, coarse_cutoff=3.5,
                 output_file="results.jsonl", csv_file=None, cluster=False, cluster_threshold=3.5):
    """
    Iterates over the trajectory and aggregates network pathways.
    Streams output as JSON Lines (JSONL) to prevent memory exhaustion,
    and optionally writes to CSV.
    """
    logger.info(f"Loading topology: {topo_file}")

    # Load universe
    if traj_file:
        u = mda.Universe(topo_file, traj_file)
    else:
        u = mda.Universe(topo_file)

    n_frames = len(u.trajectory)

    if n_frames > 1000 and stride == 1:
        logger.warning(
            "Trajectory has >1000 frames. This might consume significant time. "
            "Consider using a larger --stride parameter or clustering your trajectory."
        )

    frames_to_process = range(0, n_frames, stride)
    total_processed = len(frames_to_process)

    logger.info(f"Processing {total_processed} frames out of {n_frames} with stride {stride}.")

    # Pre-compute selections
    water_atoms = u.select_atoms(water_sel)
    root_atoms = u.select_atoms(root_sel)

    if len(root_atoms) == 0:
        logger.error(f"Error: Root selection '{root_sel}' returned no atoms. Check your resname/resid.")
        sys.exit(1)

    # Global stats
    total_length_sum = 0
    total_prob_sum = 0
    total_paths = 0
    found_large_path = False

    path_frequency = {}

    out_f = open(output_file, 'w')

    csv_f = None
    csv_writer = None
    if csv_file:
        csv_f = open(csv_file, 'w', newline='')
        csv_writer = csv.writer(csv_f)
        csv_writer.writerow(["Frame", "Root_Residue", "Path_Atom_Indices", "Path_Length", "Cumulative_Probability", "Average_OO_Distance"])

    # Metadata as the first line of JSONL
    metadata = {
        "type": "metadata",
        "n_frames_analyzed": total_processed,
        "parameters": {
            "root_sel": root_sel,
            "water_sel": water_sel,
            "stride": stride,
            "max_depth": max_depth,
            "prob_threshold": prob_threshold,
            "coarse_cutoff": coarse_cutoff
        }
    }
    out_f.write(json.dumps(metadata) + '\n')

    for ts in u.trajectory[::stride]:
        frame_idx = ts.frame
        logger.info(f"Analyzing frame {frame_idx}...")

        # 1. Build coarse graph
        g, root_indices = build_graph(u, water_atoms, root_atoms, max_distance=coarse_cutoff, max_depth=max_depth)

        # 2. Compute fine probabilities
        g = compute_edge_probabilities(g, u)

        # 3. Traverse
        paths = traverse_network(g, root_indices, max_depth=max_depth, prob_threshold=prob_threshold)

        total_paths += len(paths)

        frame_paths_data = []
        for path_indices, prob in paths:
            path_len = len(path_indices) - 1

            if path_len < min_depth:
                continue

            if path_len >= 8:
                found_large_path = True

            total_length_sum += path_len
            total_prob_sum += prob

            path_tuple = tuple(int(n) for n in path_indices)
            path_frequency.setdefault(path_tuple, set()).add(frame_idx)

            coords = []
            distances = []
            for i, node_idx in enumerate(path_indices):
                coords.append(u.atoms[node_idx].position.tolist())
                if i > 0:
                    distances.append(g[path_indices[i-1]][node_idx]['dist'])

            avg_oo = float(np.mean(distances)) if distances else 0.0

            # Also store explicit atom ids for Chimera indexing
            atom_ids = [int(u.atoms[n].id) for n in path_indices]

            frame_paths_data.append({
                "nodes": [int(n) for n in path_indices],
                "atom_ids": atom_ids,
                "coords": coords,
                "probability": float(prob),
                "length": path_len,
                "avg_oo_dist": avg_oo
            })

            if csv_writer:
                root_res = u.atoms[path_indices[0]].resname
                path_str = "-".join(str(n) for n in path_indices)

                # Sanitize to prevent CSV formula injection
                safe_root_res = sanitize_csv_field(root_res)
                safe_path_str = sanitize_csv_field(path_str)

                csv_writer.writerow([frame_idx, safe_root_res, safe_path_str, path_len, prob, avg_oo])

        # Write frame data immediately to release RAM
        out_f.write(json.dumps({"type": "frame", "frame_idx": frame_idx, "paths": frame_paths_data}) + '\n')

    # Calculate path statistics
    top_paths = []
    for p_tuple, frames_set in path_frequency.items():
        frame_count = len(frames_set)
        occupancy = float(frame_count) / total_processed if total_processed > 0 else 0.0
        top_paths.append({
            "nodes": list(p_tuple),
            "frame_count": frame_count,
            "occupancy": occupancy
        })

    # Sort descending by frame count and select top 20
    top_paths.sort(key=lambda x: x["frame_count"], reverse=True)
    top_paths = top_paths[:20]

    path_stats_obj = {
        "type": "path_statistics",
        "total_frames_analyzed": total_processed,
        "top_paths": top_paths
    }
    out_f.write(json.dumps(path_stats_obj) + '\n')

    out_f.close()
    # Check threshold notification
    if max_depth >= 8 and not found_large_path:
        notice_msg = "Notice: No water-bridges larger than 8 layers were detected in this analysis."
        logger.info(notice_msg)
        with open(output_file, 'a') as f:
            f.write(json.dumps({"type": "notice", "message": notice_msg}) + '\n')

    if csv_f:
        csv_f.close()
        logger.info(f"CSV saved to {csv_file}")

    avg_length = float(total_length_sum) / total_paths if total_paths > 0 else 0.0
    avg_prob = float(total_prob_sum) / total_paths if total_paths > 0 else 0.0

    logger.info("=== Analysis Complete ===")
    logger.info(f"Total paths found across analyzed frames: {total_paths}")
    logger.info(f"Average path length (depth): {avg_length:.2f}")
    logger.info(f"Average cumulative probability: {avg_prob:.4f}")
    logger.info(f"Results saved to {output_file}")

    if cluster:
        cluster_pathways(data_file=output_file, threshold=cluster_threshold)

    return None
