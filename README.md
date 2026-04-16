# Single-point-Water-bridge
Perform a comparative and integrative analysis for a new hybrid algorithm that adapts tunnel-search concepts from CAVER to water bridge detection.

**Objective:** Develop a standalone Python package that implements a hybrid algorithm merging CAVER's tunnel-search logic with continuous water bridge detection to explore open-ended, water-mediated hydrogen-bond networks in molecular dynamics trajectories.

**Technical Specifications & Algorithm Design:**

1.  **Framework & Data Parsing:** * Use `MDAnalysis` to parse topology and trajectory data, handle periodic boundary conditions, and perform initial coarse distance filtering between selections.
    * Use `NetworkX` to construct the adjacency matrix and execute graph traversal algorithms.

2.  **Graph Representation:** * Nodes: Individual water molecules and protein/ligand hydrogen-bond donors or acceptors.
    * Edges: Hydrogen bonds connecting these nodes.

3.  **Search Strategy:** * Implement Dijkstra's algorithm originating from a single, user-defined root coordinate (e.g., a specific catalytic residue or ligand atom). 
    * The algorithm must explore outward radially without requiring a predefined target endpoint, diverging from standard shortest-path implementations.

4.  **Edge Evaluation (Two-Stage Filter):**
    * *Coarse Filter:* Apply a strict geometric maximum heavy-atom distance cutoff (e.g., 3.5 Å) to rapidly define potential interacting pairs.
    * *Fine Evaluation (Probabilistic):* Replace binary angle/distance cutoffs with continuous mathematical switching functions. Calculate a fractional connection probability ($P_i$) for each hydrogen bond based on optimal geometries.

5.  **Pruning Rules:**
    * *Depth Limit:* Enforce a maximum path length (number of intermediate waters) to prevent infinite loops into the bulk solvent.
    * *Energy/Probability Threshold:* Terminate branch expansion if the cumulative path probability ($\prod P_i$) falls below a user-defined minimum threshold.

6.  **Scoring Function:** * Rank identified pathways using a path cost function defined as $-\ln(\prod P_{i})$. 
    * This formulation penalizes longer paths while rewarding highly stable connections, mirroring the theoretical basis of CAVER's continuous throughput calculation ($\int_{0}^{L}r(l)^{-2}dl$)[cite: 1].

**Deliverables:**
* Modular Python source code (separating parsing, graph construction, pathfinding, and scoring logic).
* A Command-Line Interface (CLI) accepting topology files, trajectory files, the root coordinate, and configurable thresholds (max depth, probability cutoff).
* Export functionality to output the ranked pathways in a format compatible with visualization tools like VMD or PyMOL (e.g., generating pseudo-atoms or connection scripts).
