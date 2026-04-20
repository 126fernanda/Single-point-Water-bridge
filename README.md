# water_bridges-nw

`water_bridges-Nw` is a high-performance Python framework for the discovery and statistical analysis of water-mediated hydrogen-bond networks (water bridges) in molecular dynamics (MD) trajectories. By merging the directional logic of tunnel-search algorithms with a continuous probabilistic model, it identifies functional "water wires" and multi-order water bridge pathways that traditional geometric cut-offs often miss.

## Features and Utilities

*   **Multi-Order Water Bridge Discovery:** Unlike tools limited to first-order bridges, this package utilizes a shell-based Breadth-First Search (BFS) to identify complex, higher-order water bridges. It originates from a single user-defined root to explore how these bridges organically extend into the bulk solvent.
*   **Continuous Probabilistic Bridge Scoring:** Replaces rigid binary cut-offs with fractional switching functions. This evaluates the stability of a water bridge based on:
    *   **Heavy-Atom Geometry:** Tunable $r_{OO}$ distance factors.
    *   **Angular Alignment:** Rewards linear donor-hydrogen-acceptor arrangements.
    *   **Chemical Specificity:** Adjusts bridge probability based on atom elements (e.g., Nitrogen vs. Oxygen vs. Sulfur).
*   **Dijkstra-Driven Path Optimization:** Employs a weighted graph traversal to determine the "Most Probable Path" for any given water bridge network. By converting probabilities into logarithmic weights, it effectively prunes low-probability noise to focus on chemically significant water-mediated interactions.
*   **Temporal Bridge Clustering:** Includes a post-processing module that uses Directed Hausdorff Distances to perform hierarchical clustering of pathways. This allows for the identification of "Collective Pathways" persistent bridge networks that appear across multiple frames, calculating their occupancy and defining representative clusters.
*   **Two-Phase Execution Architecture:**
    *   **Phase 1 (Calculate):** Processes trajectory data frame-by-frame using `MDAnalysis` and `NetworkX` to evaluate potential edges and aggregate statistical pathways. Employs a zero-memory-growth O(1) tracking architecture, streaming outputs safely via `JSON Lines (.jsonl)` and structured `.csv` to prevent RAM exhaustion.
    *   **Phase 2 (Visualize):** Parses the calculated data and seamlessly generates visualization scripts.
*   **High-Performance Visualizations:** Decoupled visualization logic exports powerful scripts targeting PyMOL (CGO cylinders), VMD (Dynamic Tcl selection indexing), and UCSF Chimera (`.py`/`.cmd`). Explicitly maps the true topological `.id` to guarantee 100% accurate visual rendering of pathways even across structurally gapped or complex topologies.

## Installation

The package requires Python 3.8 or higher. You can install `water_bridges-Nw` directly via pip.

To install it from the repository root:

```bash
pip install .
```

For development mode, run:

```bash
pip install -e .
```

## Usage

The package is run via a unified Command-Line Interface (CLI): `water_bridges_nw`.

### 1. Calculation Phase
Analyze your MD trajectory to discover water-mediated networks. By default, the output is saved to `results.json`.

```bash
water_bridges_nw calculate \
  --topo my_topology.pdb \
  --traj my_trajectory.xtc \
  --root "resname LIG and name O1" \
  --stride 10 \
  --max_depth 10 \
  --min_depth 3 \
  --prob_threshold 0.001 \
  --output custom_results.jsonl \
  --csv custom_summary.csv
```

*Options:*
*   `--topo`: Topology file (.pdb, .tpr, etc.)
*   `--traj`: Trajectory file (.xtc, .dcd). If omitted, evaluates only the topology.
*   `--root`: Standard MDAnalysis atom selection string defining the starting root coordinate. Supports `resname`, `resid`, and `name` strings.
*   `--water`: Selection string for solvent (default: `"resname SOL or resname WAT or resname HOH"`).
*   `--stride`: Frame stride to process. *(Note: A warning is issued if the evaluated frames exceed 1000)*
*   `--max_depth`: Maximum chain length of sequential waters (default: `10`).
*   `--min_depth`: Minimum chain length filter. Shorter paths will be discarded (default: `1`).
*   `--prob_threshold`: Cumulative probabilistic threshold; paths falling below this probability are pruned.
*   `--output`: Custom filename for the detailed `JSON Lines` output (default: `results.jsonl`).
*   `--csv`: Optional custom filename to generate a structured, human-readable summary of the detected paths.

### 2. Visualization Phase
Generate rendering scripts from your calculation output.

**Generate a density map (all frames overlaid) for PyMOL:**
```bash
water_bridges_nw visualize \
  --data custom_results.jsonl \
  --format pymol \
  --mode density \
  --output my_density_network.py
```

**Generate a pathway view for a specific frame in VMD:**
```bash
water_bridges_nw visualize \
  --data custom_results.jsonl \
  --format vmd \
  --mode frame \
  --frame 10 \
  --output frame_10_network.tcl
```

**Generate a pathway view for a specific frame in UCSF Chimera:**
```bash
water_bridges_nw visualize \
  --data custom_results.jsonl \
  --format chimera \
  --mode frame \
  --frame 10 \
  --output frame_10_network.py
```

*Options:*
*   `--data`: The `.jsonl` output generated by the `calculate` subcommand.
*   `--format`: Target visualization software (`vmd`, `pymol`, or `chimera`).
*   `--mode`: Either `density` (overlays pathways across all analyzed frames) or `frame` (extracts a single frame).
*   `--frame`: Index of the frame to visualize (required if `--mode frame` is used).
*   `--output`: Custom filename/prefix for the exported script to prevent overwriting outputs (default: `pathways_viz`).
