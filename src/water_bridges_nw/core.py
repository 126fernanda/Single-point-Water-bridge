import numpy as np
import MDAnalysis as mda
import networkx as nx
from MDAnalysis.lib.distances import capped_distance, calc_bonds
from .math_utils import calculate_hbond_probability

def build_graph(u, water_atoms, root_atoms, max_distance=3.5, max_depth=5):
    """
    Builds a NetworkX graph representing potential hydrogen bonds
    using a shell-based iterative expansion to avoid global N^2 distance matrices.
    """
    g = nx.Graph()

    for a in root_atoms:
        g.add_node(a.index, resname=a.resname, name=a.name, pos=a.position)

    current_shell_indices = list(root_atoms.indices)
    visited_indices = set(current_shell_indices)

    box = u.dimensions

    for depth in range(max_depth):
        if not current_shell_indices:
            break

        current_shell = u.atoms[current_shell_indices]

        # capped_distance between current shell and all waters
        pairs, distances = capped_distance(
            current_shell.positions,
            water_atoms.positions,
            max_cutoff=max_distance,
            box=box,
            return_distances=True
        )

        next_shell_indices = set()

        for (idx_current, idx_water), dist in zip(pairs, distances):
            u_node = current_shell.indices[idx_current]
            v_node = water_atoms.indices[idx_water]

            if v_node not in g:
                a2 = u.atoms[v_node]
                g.add_node(v_node, resname=a2.resname, name=a2.name, pos=a2.position)

            g.add_edge(u_node, v_node, dist=dist)

            if v_node not in visited_indices:
                next_shell_indices.add(v_node)

        # Connect nodes within the new shell
        if next_shell_indices:
            next_shell_list = list(next_shell_indices)
            next_shell_ag = u.atoms[next_shell_list]
            intra_pairs, intra_dist = capped_distance(
                next_shell_ag.positions,
                next_shell_ag.positions,
                max_cutoff=max_distance,
                box=box,
                return_distances=True
            )
            for (i1, i2), d in zip(intra_pairs, intra_dist):
                if i1 < i2:
                    n1 = next_shell_ag.indices[i1]
                    n2 = next_shell_ag.indices[i2]
                    g.add_edge(n1, n2, dist=d)

        visited_indices.update(next_shell_indices)
        current_shell_indices = list(next_shell_indices)

    return g, root_atoms.indices

def compute_edge_probabilities(g, u):
    """
    Iterates over edges in the graph and computes the continuous probability
    by dynamically finding attached hydrogens.
    """
    edges_to_remove = []
    h_cache = {}

    def get_hydrogens(atom):
        if atom.index in h_cache:
            return h_cache[atom.index]
        hs = [a for a in atom.residue.atoms if a.name.startswith('H') or a.type == 'H']
        h_cache[atom.index] = hs
        return hs

    for u_node, v_node, data in g.edges(data=True):
        a1 = u.atoms[u_node]
        a2 = u.atoms[v_node]

        mod_rOO = data['dist']

        hs1 = get_hydrogens(a1)
        hs2 = get_hydrogens(a2)
        all_hs = hs1 + hs2

        if not all_hs:
            # Fallback
            mod_rOiH = 1.0
            mod_rOjH = mod_rOO - 1.0
        else:
            min_diff = float('inf')
            best_rOiH = 1.0
            best_rOjH = mod_rOO - 1.0

            p1 = a1.position
            p2 = a2.position

            for h in all_hs:
                ph = h.position
                if u.dimensions is not None:
                    d1 = calc_bonds(p1, ph, box=u.dimensions)
                    d2 = calc_bonds(p2, ph, box=u.dimensions)
                else:
                    d1 = np.linalg.norm(p1 - ph)
                    d2 = np.linalg.norm(p2 - ph)

                if (d1 + d2) < min_diff:
                    min_diff = d1 + d2
                    best_rOiH = d1
                    best_rOjH = d2

            mod_rOiH = best_rOiH
            mod_rOjH = best_rOjH

        prob = calculate_hbond_probability(mod_rOO, mod_rOiH, mod_rOjH)

        if prob <= 0:
            edges_to_remove.append((u_node, v_node))
        else:
            score = -np.log(prob) if prob > 0 else float('inf')
            g[u_node][v_node]['prob'] = prob
            g[u_node][v_node]['weight'] = score

    g.remove_edges_from(edges_to_remove)
    return g

def traverse_network(g, root_indices, max_depth=5, prob_threshold=1e-3):
    """
    Performs Dijkstra-based search to find the optimal path to each reachable node.
    """
    import heapq

    dist = {}
    best_path = {}

    pq = []

    for root in root_indices:
        if root in g:
            heapq.heappush(pq, (0.0, root, 0, [root]))
            dist[root] = 0.0
            best_path[root] = [root]

    while pq:
        curr_weight, u_node, depth, path = heapq.heappop(pq)

        if curr_weight > dist.get(u_node, float('inf')):
            continue

        if depth >= max_depth:
            continue

        for v_node in g.neighbors(u_node):
            if v_node in path:
                continue

            edge_weight = g[u_node][v_node]['weight']
            next_weight = curr_weight + edge_weight
            next_prob = np.exp(-next_weight)

            if next_prob >= prob_threshold:
                if next_weight < dist.get(v_node, float('inf')):
                    dist[v_node] = next_weight
                    next_path = path + [v_node]
                    best_path[v_node] = next_path
                    heapq.heappush(pq, (next_weight, v_node, depth + 1, next_path))

    final_paths = []
    for node, path in best_path.items():
        if len(path) > 1:
            prob = np.exp(-dist[node])
            if prob >= prob_threshold:
                final_paths.append((path, float(prob)))

    return final_paths
