import json
import logging
import csv
from collections import defaultdict
import MDAnalysis as mda
import numpy as np

from .core import build_graph, compute_edge_probabilities, traverse_network

logger = logging.getLogger(__name__)

def run_analysis(topo_file, traj_file, root_sel, water_sel="resname SOL or resname WAT or resname HOH",
                 stride=1, max_depth=5, prob_threshold=1e-3, coarse_cutoff=3.5,
                 output_file="results.jsonl", csv_file=None):
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

    # Global stats
    total_length_sum = 0
    total_prob_sum = 0
    total_paths = 0

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
            total_length_sum += path_len
            total_prob_sum += prob

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
                csv_writer.writerow([frame_idx, root_res, path_str, path_len, prob, avg_oo])

        # Write frame data immediately to release RAM
        out_f.write(json.dumps({"type": "frame", "frame_idx": frame_idx, "paths": frame_paths_data}) + '\n')

    out_f.close()
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

    # We could return a dictionary for backwards compatibility but we rely on files now.
    return None
