import numpy as np
import MDAnalysis as mda
import networkx as nx
from MDAnalysis.lib.distances import capped_distance, calc_bonds, distance_array
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

        candidate_hs = [a for a in atom.residue.atoms if a.name.startswith('H') or a.type == 'H']
        if not candidate_hs:
            h_cache[atom.index] = []
            return []

        hs_positions = np.array([h.position for h in candidate_hs])
        dists = calc_bonds(atom.position, hs_positions, box=u.dimensions)

        bonded_hs = [candidate_hs[i] for i, d in enumerate(dists) if d <= 1.2]

        h_cache[atom.index] = bonded_hs
        return bonded_hs

    for u_node, v_node, data in g.edges(data=True):
        a1 = u.atoms[u_node]
        a2 = u.atoms[v_node]

        mod_rOO = data['dist']

        hs1 = get_hydrogens(a1)
        hs2 = get_hydrogens(a2)
        all_hs = hs1 + hs2

        if not all_hs:
            prob = 0.0
        else:
            p_hs = np.array([h.position for h in all_hs])
            d1_array = distance_array(np.array([a1.position]), p_hs, box=u.dimensions)[0]
            d2_array = distance_array(np.array([a2.position]), p_hs, box=u.dimensions)[0]

            # Apply strict geometric filters:
            # - Distance cutoff: At least one of the OH distances (donor-hydrogen) must be reasonably short (e.g., covalent bond ~ 1.0A).
            # - Acceptor-Hydrogen distance <= 3.0 A
            # - Donor-Hydrogen-Acceptor Angle >= 120 degrees
            # Since a1 and a2 are heavy atoms (O, N, etc.), one acts as donor, one as acceptor.

            from MDAnalysis.lib.distances import calc_angles

            best_prob = 0.0

            for idx, h_pos in enumerate(p_hs):
                # Is a1 the donor or a2?
                # The covalent bond is typically < 1.2 A
                is_a1_donor = d1_array[idx] < 1.2
                is_a2_donor = d2_array[idx] < 1.2

                if not (is_a1_donor or is_a2_donor):
                    continue

                donor_pos = a1.position if is_a1_donor else a2.position
                acceptor_pos = a2.position if is_a1_donor else a1.position

                # Acceptor-Hydrogen distance
                dist_HA = d2_array[idx] if is_a1_donor else d1_array[idx]
                if dist_HA > 3.0:
                    continue

                # Angle: Donor - Hydrogen - Acceptor
                # calc_angles takes arrays of coords: (pos1, pos2, pos3) for angle 1-2-3
                # We need angle donor - hydrogen - acceptor
                # Note: MDAnalysis calc_angles expects shapes (N, 3)
                angle = calc_angles(
                    np.array([donor_pos]),
                    np.array([h_pos]),
                    np.array([acceptor_pos]),
                    box=u.dimensions
                )[0]

                # Convert radians to degrees (calc_angles returns radians)
                angle_deg = np.degrees(angle)

                if angle_deg >= 120.0:
                    # Determine appropriate r0_oo based on heavy atom elements
                    try:
                        e1 = a1.element
                    except AttributeError:
                        e1 = a1.name[0]
                    try:
                        e2 = a2.element
                    except AttributeError:
                        e2 = a2.name[0]

                    e1 = e1.upper()
                    e2 = e2.upper()

                    if 'S' in (e1, e2):
                        r0_oo = 3.3
                    elif (e1, e2) in (('N', 'N'),):
                        r0_oo = 3.0
                    elif (e1, e2) in (('N', 'O'), ('O', 'N')):
                        r0_oo = 2.9
                    else: # O-O and defaults
                        r0_oo = 2.7

                    # Valid H-bond geometry. Calculate continuous probability.
                    mod_rOiH = d1_array[idx]
                    mod_rOjH = d2_array[idx]
                    p = calculate_hbond_probability(mod_rOO, mod_rOiH, mod_rOjH, r0_oo=r0_oo)
                    if p > best_prob:
                        best_prob = p

            prob = best_prob

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
