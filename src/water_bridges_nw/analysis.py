import json
import logging
from collections import defaultdict
import MDAnalysis as mda
import numpy as np

from .core import build_graph, compute_edge_probabilities, traverse_network

logger = logging.getLogger(__name__)

def run_analysis(topo_file, traj_file, root_sel, water_sel="resname SOL or resname WAT or resname HOH",
                 stride=1, max_depth=5, prob_threshold=1e-3, coarse_cutoff=3.5,
                 output_file="results.json"):
    """
    Iterates over the trajectory and aggregates network pathways.
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
            "Trajectory has >1000 frames. This might consume significant time and memory. "
            "Consider using a larger --stride parameter or clustering your trajectory."
        )

    frames_to_process = range(0, n_frames, stride)
    total_processed = len(frames_to_process)

    logger.info(f"Processing {total_processed} frames out of {n_frames} with stride {stride}.")

    # We want to aggregate statistics about paths
    # A path is a sequence of atom types or relative positions, but exact indices change if waters diffuse.
    # For meaningful aggregation, we can store properties of paths found per frame.
    # For visualization, we will store the raw cylinder coordinates or node coordinates for frames.

    # Pre-compute selections to avoid string parsing every frame
    water_atoms = u.select_atoms(water_sel)
    root_atoms = u.select_atoms(root_sel)

    # We will store results per frame for visualization
    results_by_frame = {}

    # Global stats
    all_path_lengths = []
    all_path_probs = []
    total_paths = 0

    for ts in u.trajectory[::stride]:
        frame_idx = ts.frame
        logger.info(f"Analyzing frame {frame_idx}...")

        # 1. Build coarse graph
        g, root_indices = build_graph(u, water_atoms, root_atoms, max_distance=coarse_cutoff)

        # 2. Compute fine probabilities
        g = compute_edge_probabilities(g, u, water_sel, root_sel)

        # 3. Traverse
        paths = traverse_network(g, root_indices, max_depth=max_depth, prob_threshold=prob_threshold)

        total_paths += len(paths)

        frame_paths_data = []
        for path_indices, prob in paths:
            all_path_lengths.append(len(path_indices) - 1)
            all_path_probs.append(prob)

            # Store coordinates of the nodes in the path for visualization
            coords = []
            for node_idx in path_indices:
                coords.append(u.atoms[node_idx].position.tolist())

            frame_paths_data.append({
                "nodes": [int(n) for n in path_indices],
                "coords": coords,
                "probability": float(prob),
                "length": len(path_indices) - 1
            })

        results_by_frame[frame_idx] = frame_paths_data

    # Aggregate statistics
    avg_length = np.mean(all_path_lengths) if all_path_lengths else 0.0
    avg_prob = np.mean(all_path_probs) if all_path_probs else 0.0

    logger.info("=== Analysis Complete ===")
    logger.info(f"Total paths found across analyzed frames: {total_paths}")
    logger.info(f"Average path length (depth): {avg_length:.2f}")
    logger.info(f"Average cumulative probability: {avg_prob:.4f}")

    summary = {
        "metadata": {
            "n_frames_analyzed": total_processed,
            "total_paths": total_paths,
            "avg_length": float(avg_length),
            "avg_probability": float(avg_prob),
            "parameters": {
                "root_sel": root_sel,
                "water_sel": water_sel,
                "stride": stride,
                "max_depth": max_depth,
                "prob_threshold": prob_threshold,
                "coarse_cutoff": coarse_cutoff
            }
        },
        "frames": results_by_frame
    }

    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Results saved to {output_file}")
    return summary
