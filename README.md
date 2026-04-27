# water_bridges-nw

`water_bridges-nw` is a Python framework for the discovery and statistical
analysis of water-mediated hydrogen-bond networks (water bridges) in molecular
dynamics (MD) trajectories. It combines a continuous probabilistic H-bond
scoring model with graph-based path enumeration to identify multi-order water
bridge pathways and their temporal occupancy across simulation frames.

## Features

- **Multi-order water bridge discovery:** A depth-bounded graph traversal
  originating from a user-defined root residue identifies water bridges
  up to `max_depth` consecutive waters into the solvent.

- **Continuous probabilistic scoring:** Replaces binary geometric cutoffs
  with fractional switching functions that evaluate:
  - Heavy-atom O···O distance
  - Angular donor-hydrogen-acceptor alignment via a triangle-inequality proxy
  - Steric repulsion for clash geometries (r_OO < 2.40 Å)
  - Chemical specificity: per-element r0 reference distances for N, O, S

- **Weighted path enumeration:** Edge probabilities are converted to
  logarithmic weights so that path traversal naturally ranks higher-probability
  chains above weaker ones within the depth limit.

- **Temporal occupancy statistics:** Tracks which specific water-molecule
  chains appear in each frame and reports the top 20 paths ranked by
  frame-count occupancy.

- **Pathway clustering:** Optional post-processing groups spatially similar
  paths across frames using directed Hausdorff distance and average-link
  hierarchical clustering, producing a `clustered_pathways.json` file with
  cluster size, occupancy, average probability, and a representative
  medoid geometry.

- **Two-phase execution:**
  - **Phase 1 — `calculate`:** Processes trajectory frames via MDAnalysis
    and NetworkX. Streams output as JSON Lines (`.jsonl`) and optional `.csv`
    to avoid RAM exhaustion on long trajectories.
  - **Phase 2 — `visualize`:** Parses `.jsonl` output and generates
    ready-to-run scripts for PyMOL (CGO cylinders), VMD (Tcl atom
    index selection), and UCSF Chimera (`.py` / `.bild`).

## Installation

Requires Python 3.8 or higher.

```bash
# Standard install from repository root
pip install .

# Development (editable) install
pip install -e .
```

## Usage

### 1. Calculate

```bash
water_bridges_nw calculate \
  --topo  my_topology.pdb \
  --traj  my_trajectory.xtc \
  --root  "resname LIG and name O1" \
  --stride 10 \
  --max_depth 10 \
  --min_depth 3 \
  --coarse_cutoff 4.5 \
  --output results.jsonl \
  --csv   summary.csv \
  --cluster \
  --cluster_threshold 3.5
```

| Option | Default | Description |
|---|---|---|
| `--topo` | required | Topology file (.pdb, .tpr, …) |
| `--traj` | — | Trajectory file (.xtc, .dcd, …). Omit to evaluate topology only. |
| `--root` | required | MDAnalysis selection string for the root atom(s) (e.g. `"resname LIG"`). |
| `--water` | `"resname SOL or resname WAT or resname HOH"` | Solvent selection string. |
| `--stride` | `1` | Process every Nth frame. A warning is issued when `stride=1` and the trajectory exceeds 1000 frames. |
| `--max_depth` | `10` | Maximum number of sequential water molecules in a path. |
| `--min_depth` | `1` | Minimum path length; shorter paths are discarded from output. |
| `--coarse_cutoff` | `4.5` | Distance cutoff in Å for initial neighbour graph construction. |
| `--output` | `results.jsonl` | Output file for full frame-by-frame path data (JSON Lines format). |
| `--csv` | — | Optional human-readable CSV summary of detected paths. |
| `--cluster` | off | Enable post-analysis spatial clustering of paths. Writes `clustered_pathways.json`. |
| `--cluster_threshold` | `3.5` | Hausdorff distance threshold in Å for clustering. |


### 2. Visualize

```bash
# Density overlay across all frames — PyMOL
water_bridges_nw visualize \
  --data results.jsonl --format pymol --mode density \
  --output network_density.py

# Single-frame selection — VMD
water_bridges_nw visualize \
  --data results.jsonl --format vmd --mode frame --frame 10 \
  --output frame_10.tcl

# Single-frame selection — UCSF Chimera
water_bridges_nw visualize \
  --data results.jsonl --format chimera --mode frame --frame 10 \
  --output frame_10.py
```

| Option | Description |
|---|---|
| `--data` | `.jsonl` file produced by `calculate`. |
| `--format` | `vmd`, `pymol`, or `chimera`. |
| `--mode` | `density` (all frames overlaid) or `frame` (single frame). |
| `--frame` | Frame index to visualize; required when `--mode frame`. |
| `--output` | Output script filename or prefix (default: `pathways_viz`). |

## Solvent naming conventions

The default water selection covers GROMACS (`SOL`), AMBER (`WAT`), and
CHARMM/PDB (`HOH`) residue names. For non-standard solvents or co-solvents
acting as bridge donors, pass a custom `--water` string, for example:

```bash
--water "resname SOL or resname GOL"
```

For united-atom force fields (GROMOS, OPLS-UA) or coarse-grained models,
consult the force-field documentation for the correct oxygen atom names
and verify that hydrogens are present in the topology.

## Limitations

Understanding what this tool does and does not model helps you interpret
its output correctly.

### What the tool is designed for

- Identifying **water bridge networks** connecting a root residue to bulk
  solvent across one or more water molecules in standard all-atom MD.
- Ranking bridges by **temporal occupancy** (fraction of frames in which
  a given atom-index chain appears).
- Generating **spatial cluster representatives** for bridges that recur
  across many frames.
- Producing **density maps** of where water-mediated interactions occur
  near a binding site or protein surface.

### What the tool does not model

**Periodic boundary conditions in path coordinates.** The neighbour graph
is built with distance-based selection (which MDAnalysis applies with PBC
awareness), but the 3D coordinates stored for each path node are raw
wrapped positions. Paths that cross a periodic boundary will display as
broken or elongated segments in visualization software. This affects
trajectories where the root residue and solvent are in different periodic
images.

**Water molecule exchange along a pathway.** Occupancy is tracked by
exact atom-index tuples. If the same physical channel is traversed by
different water molecules at different times — common in trajectories
longer than a few nanoseconds — each unique permutation is counted as a
separate path. Temporal occupancy will be underestimated for highly
dynamic pathways. The spatial clustering step partially compensates
for this.

**Proton-conducting (Grotthuss) water wires.** The algorithm is
undirected: it detects chains of mutually compatible H-bond geometries
but does not verify that water dipoles are aligned head-to-tail, which
is the physical requirement for proton transport. An antiparallel pair
that blocks conductance is geometrically indistinguishable from a
conducting pair in the current model. If you are studying aquaporins,
gramicidin channels, or other systems where proton-wire directionality
matters, the occupancy output should be interpreted as a measure of
**structural presence**, not **transport competence**.

**Transmembrane pathways crossing the periodic boundary.** Wires
spanning the full membrane thickness cross periodic images. Because
path coordinates are not minimum-image corrected, these wires will
not be reliably detected or visualised.

**Consecutive-frame persistence.** A path that appears in 100 frames
scattered across a 10 000-frame trajectory receives the same
occupancy score (0.01) as a path present in 100 consecutive frames.
Only the latter constitutes a genuinely persistent structural feature.
If persistence matters for your analysis, filter the `top_paths` output
manually or apply the `--cluster` option, which groups paths by spatial
similarity and reports per-cluster occupancy.

**Very large path counts.** The clustering step computes a full
pairwise Hausdorff distance matrix, which scales as O(N²) in the
number of paths. For long trajectories with a high branching factor,
the clustering step can become slow or memory-intensive. Use
`--stride` to reduce frame count or increase `--min_depth` to
reduce the number of short paths before enabling `--cluster`.
