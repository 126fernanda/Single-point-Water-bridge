import numpy as np
import MDAnalysis as mda
import networkx as nx
from MDAnalysis.lib.distances import capped_distance
from .math_utils import calculate_hbond_probability

def build_graph(u, water_atoms, root_atoms, max_distance=3.5):
    """
    Builds a NetworkX graph representing potential hydrogen bonds.
    """
    g = nx.Graph()

    # We now expect water_atoms and root_atoms to be passed directly
    all_atoms = water_atoms + root_atoms

    # Use capped_distance to find all pairs within max_distance
    # This is our "coarse filter"
    box = u.dimensions
    pairs, distances = capped_distance(
        all_atoms.positions,
        all_atoms.positions,
        max_cutoff=max_distance,
        box=box,
        return_distances=True
    )

    # We map back from pair indices to actual atom indices
    atom_indices = all_atoms.indices

    for (idx1, idx2), dist in zip(pairs, distances):
        if idx1 >= idx2: # undirected graph, ignore self-loops and duplicates
            continue

        atom1 = all_atoms[idx1]
        atom2 = all_atoms[idx2]

        # Only add edge if it's between a root atom and a water, or water and water
        # (Assuming root_atoms and water_atoms don't overlap, and we only want valid h-bond pairs)
        # In a more advanced version, we should identify O/N/S/F atoms and hydrogens explicitly

        # We need coordinates to compute probabilities.
        # However, math_utils.calculate_hbond_probability requires distances like mod_rOO, mod_rOiH, mod_rOjH
        # This implies we need to know where the hydrogens are.
        # For simplicity in this first draft, if we don't have explicit H, we might approximate
        # or we need the user to provide an explicit H topology.
        # Let's assume a simplified proxy for the fine filter if exact H positions aren't parsed perfectly yet,
        # but the prompt implies we evaluate `calculate_hbond_probability` which requires these distances.

        # For now, let's just add the coarse edge. The actual weights will be calculated in a refined step.
        g.add_node(atom1.index, resname=atom1.resname, name=atom1.name, pos=atom1.position)
        g.add_node(atom2.index, resname=atom2.resname, name=atom2.name, pos=atom2.position)
        g.add_edge(atom1.index, atom2.index, dist=dist)

    return g, root_atoms.indices

def compute_edge_probabilities(g, u, water_sel, root_sel):
    """
    Iterates over edges in the graph and computes the continuous probability.
    """
    edges_to_remove = []

    # To compute proper rOiH, rOjH, we ideally need to find the hydrogens attached to the O atoms.
    # In MDAnalysis, we can find bonded atoms if topology has bonds.
    # If not, distance-based H-finding might be needed.

    for u_node, v_node, data in g.edges(data=True):
        a1 = u.atoms[u_node]
        a2 = u.atoms[v_node]

        # This requires detailed knowledge of the atoms.
        # Let's mock the mod_rOO, mod_rOiH, mod_rOjH calculation for now,
        # assuming a1 and a2 are the heavy atoms.
        # If we just have O-O distance:
        mod_rOO = data['dist']

        # MOCKING:
        # In a real scenario, we'd find the H's bonded to a1 and a2.
        # We will assume optimal H placement for this skeleton code to compile and run:
        # H is typically ~1.0A from O. Ideal h-bond is linear, so rOH + rH...O = rOO
        mod_rOiH = 1.0 # mock
        mod_rOjH = mod_rOO - 1.0 # mock

        prob = calculate_hbond_probability(mod_rOO, mod_rOiH, mod_rOjH)

        if prob <= 0:
            edges_to_remove.append((u_node, v_node))
        else:
            # Score = -ln(P)
            score = -np.log(prob) if prob > 0 else float('inf')
            g[u_node][v_node]['prob'] = prob
            g[u_node][v_node]['weight'] = score # Using score as weight for shortest path/Dijkstra

    g.remove_edges_from(edges_to_remove)
    return g

def traverse_network(g, root_indices, max_depth=5, prob_threshold=1e-3):
    """
    Performs open-ended BFS/Dijkstra from root nodes.
    Returns discovered paths.
    """
    # A path is a list of node indices.
    # We want to explore outward up to max_depth.

    from collections import deque
    valid_paths = []

    # For each root node, start an exploration
    for start_node in root_indices:
        if start_node not in g:
            continue

        # We can use a queue for BFS: (current_path, cumulative_prob)
        queue = deque([([start_node], 1.0)])

        while queue:
            current_path, current_prob = queue.popleft()

            # If path reaches max_depth, stop expanding
            if len(current_path) - 1 >= max_depth:
                valid_paths.append((current_path, current_prob))
                continue

            last_node = current_path[-1]
            expanded = False

            for neighbor in g.neighbors(last_node):
                if neighbor not in current_path: # Avoid loops in the current path
                    edge_prob = g[last_node][neighbor]['prob']
                    new_prob = current_prob * edge_prob

                    if new_prob >= prob_threshold:
                        queue.append((current_path + [neighbor], new_prob))
                        expanded = True

            # If we didn't expand further but it's a valid path > 1 node, we can keep it as a terminal path
            if not expanded and len(current_path) > 1:
                valid_paths.append((current_path, current_prob))

    return valid_paths
