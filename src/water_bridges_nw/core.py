import re
import logging
import numpy as np
import networkx as nx
from MDAnalysis.lib.distances import capped_distance, distance_array
from .math_utils import calculate_hbond_probability, switching_function

logger = logging.getLogger(__name__)

_NAME_TO_ELEMENT = {
    'OW': 'O', 'O1': 'O', 'O2': 'O', 'OD1': 'O', 'OD2': 'O', 'OE1': 'O', 'OE2': 'O', 'OG': 'O', 'OG1': 'O', 'OH': 'O',
    'NZ': 'N', 'ND1': 'N', 'ND2': 'N', 'NE': 'N', 'NE1': 'N', 'NE2': 'N', 'NH1': 'N', 'NH2': 'N',
    'SG': 'S', 'SD': 'S'
}

_warned_united_atom = set()

def _is_hydrogen(a):
    """
    Tiered fallback to correctly identify if an atom is hydrogen.
    1. Element (most reliable)
    2. Mass (force-field agnostic)
    3. Name/Type (last resort)
    """
    if getattr(a, 'element', '') and a.element.strip().upper() == 'H':
        return True

    try:
        if 0.5 < a.mass < 2.5:
            return True
    except Exception:
        pass

    return bool(re.search(r'(?i)\bh', a.name)) or getattr(a, 'type', '') == 'H'

def _get_element(atom):
    """
    Robustly resolves the element of an atom, preventing misclassification.
    Priority 1: atom.element (MDAnalysis standard).
    Priority 2: Exact match in _NAME_TO_ELEMENT dictionary.
    Priority 3: Strip leading digits from atom.name and take the leading alphabetic substring.
    """
    valid_elements = {"O", "N", "S", "F", "CL", "BR"}
    try:
        if atom.element:
            e = atom.element.strip().upper()
            if e in valid_elements:
                return e
    except AttributeError:
        pass

    # Fallback 1: Lookup exact names for common topologies (e.g., OW -> O, NZ -> N)
    atom_name_upper = atom.name.strip().upper()
    if atom_name_upper in _NAME_TO_ELEMENT:
        return _NAME_TO_ELEMENT[atom_name_upper]

    # Fallback 2: Strip leading digits and extract the leading alphabetic substring
    match = re.search(r'^[0-9]*([A-Za-z]+)', atom.name)
    if match:
        e = match.group(1).upper()
        if e in valid_elements:
            return e

    return "UNKNOWN"

def build_graph(u, water_atoms, root_atoms, max_distance=4.5, max_depth=5):
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

        candidate_hs = [a for a in atom.residue.atoms if _is_hydrogen(a)]
        if not candidate_hs:
            h_cache[atom.index] = []
            return []

        if candidate_hs:
            hs_positions = [h.position for h in candidate_hs]
            h_cache[atom.index] = hs_positions
            return hs_positions

        # No explicit hydrogens found - United-Atom geometry projection
        try:
            bonded_heavy_atoms = atom.bonded_atoms
        except MDAnalysis.exceptions.NoDataError:
            bonded_heavy_atoms = []

        if not bonded_heavy_atoms:
            elem = _get_element(atom)
            if elem in {"O", "N", "S"} and atom.index not in _warned_united_atom:
                logger.warning(f"Heavy atom {elem} (index {atom.index}) has no bonded hydrogens and no neighbors. Virtual projection failed.")
                _warned_united_atom.add(atom.index)
            h_cache[atom.index] = []
            return []

        # Hybridization-aware United-Atom geometry projection
        neighbor_pos = np.array([neighbor.position for neighbor in bonded_heavy_atoms])

        if len(bonded_heavy_atoms) == 1:
            # sp2 approximation (e.g., carbonyl oxygen)
            v1 = atom.position - neighbor_pos[0]
            v1 /= np.linalg.norm(v1)

        if norm < 1e-6: # To avoid division by zero if positions exactly overlap
            h_cache[atom.index] = []
            return []
        else:
            # Scale to 1.0 A
            virtual_h = atom.position + (vec / norm) * 1.0
            # Generate an arbitrary orthogonal vector
            arb = np.array([1.0, 0.0, 0.0]) if abs(v1[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
            perp = np.cross(v1, arb)
            perp /= np.linalg.norm(perp)

            # Create two lone pair vectors at 120 degrees from the incoming bond
            cos_120, sin_120 = -0.5, 0.866
            lp1 = atom.position + (v1 * cos_120 + perp * sin_120) * 1.0
            lp2 = atom.position + (v1 * cos_120 - perp * sin_120) * 1.0

            h_cache[atom.index] = [lp1, lp2]
            return [lp1, lp2]

        elif len(bonded_heavy_atoms) == 2:
            # sp3 approximation (e.g., ether oxygen, thioether sulfur)
            v1 = neighbor_pos[0] - atom.position
            v2 = neighbor_pos[1] - atom.position
            v1 /= np.linalg.norm(v1)
            v2 /= np.linalg.norm(v2)

            # Normal to the plane defined by Atom and its 2 neighbors
            n = np.cross(v1, v2)
            n /= np.linalg.norm(n)

            # Vector bisecting the angle
            bisector = v1 + v2
            bisector /= np.linalg.norm(bisector)


            # Normal to the plane defined by Atom and its 2 neighbors
            n = np.cross(v1, v2)
            n /= np.linalg.norm(n)

            # Vector bisecting the angle
            bisector = v1 + v2
            bisector /= np.linalg.norm(bisector)

            # Project lone pairs outwards (tetrahedral geometry)
            cos_tilt, sin_tilt = -0.577, 0.816 # approximate projection relative to bisector
            lp1 = atom.position + (-bisector * cos_tilt + n * sin_tilt) * 1.0
            lp2 = atom.position + (-bisector * cos_tilt - n * sin_tilt) * 1.0

            h_cache[atom.index] = [lp1, lp2]
            return [lp1, lp2]


            # Normal to the plane defined by Atom and its 2 neighbors
            n = np.cross(v1, v2)
            n /= np.linalg.norm(n)

            # Vector bisecting the angle
            bisector = v1 + v2
            bisector /= np.linalg.norm(bisector)

            # Project lone pairs outwards (tetrahedral geometry)
            cos_tilt, sin_tilt = -0.577, 0.816 # approximate projection relative to bisector
            lp1 = atom.position + (-bisector * cos_tilt + n * sin_tilt) * 1.0
            lp2 = atom.position + (-bisector * cos_tilt - n * sin_tilt) * 1.0

            h_cache[atom.index] = [lp1, lp2]
            return [lp1, lp2]

        else:
            h_cache[atom.index] = []
            return []

    for u_node, v_node, data in g.edges(data=True):
        a1 = u.atoms[u_node]
        a2 = u.atoms[v_node]

        e1 = _get_element(a1)
        e2 = _get_element(a2)

        if e1 == "UNKNOWN" or e2 == "UNKNOWN":
            edges_to_remove.append((u_node, v_node))
            continue

        mod_rOO = data['dist']

        hs1 = get_hydrogens(a1)
        hs2 = get_hydrogens(a2)
        all_hs = hs1 + hs2

        # Determine appropriate r0_oo based on heavy atom elements
        if 'S' in (e1, e2):
            r0_oo_fixed = 3.3
            r0_threshold_fixed = 0.8
        elif (e1, e2) in (('N', 'N'),):
            r0_oo_fixed = 3.0
            r0_threshold_fixed = 0.6
        elif (e1, e2) in (('N', 'O'), ('O', 'N')):
            r0_oo_fixed = 2.9
            r0_threshold_fixed = 0.55
        else: # O-O and defaults
            r0_oo_fixed = 2.80
            r0_threshold_fixed = 0.45

        if not hs1 or not hs2:
            # United-atom fallback: either side is missing hydrogens, ignore angle
            prob = calculate_hbond_probability(
            best_prob = calculate_hbond_probability(
                mod_rOO, None, None,
                r0_oo=r0_oo_fixed,
                r0_threshold=r0_threshold_fixed,
                ignore_angle=True
            )
        else:
            p_hs = np.array(all_hs)
            d1_array = distance_array(np.array([a1.position]), p_hs, box=u.dimensions)[0]
            d2_array = distance_array(np.array([a2.position]), p_hs, box=u.dimensions)[0]

            best_prob = 0.0

            for idx, h_pos in enumerate(p_hs):
                dist_DH = min(d1_array[idx], d2_array[idx])
                dist_HA = max(d1_array[idx], d2_array[idx])

                p_base = calculate_hbond_probability(
                    mod_rOO, dist_DH, dist_HA,
                    r0_oo=r0_oo_fixed,
                    r0_threshold=r0_threshold_fixed
                )
                p_ha = switching_function(dist_HA, threshold=2.5, power_num=6, power_den=12)
                p_covalent = switching_function(dist_DH, threshold=1.1, power_num=6, power_den=12)

                p_i = p_base * p_ha * p_covalent
                best_prob = 1.0 - (1.0 - best_prob) * (1.0 - p_i)

        if best_prob <= 0:
            edges_to_remove.append((u_node, v_node))
        else:
            score = -np.log(best_prob) if best_prob > 0 else float('inf')
            g[u_node][v_node]['prob'] = best_prob
            g[u_node][v_node]['weight'] = score

    g.remove_edges_from(edges_to_remove)

    return g

def traverse_network(g, root_indices, max_depth=5, prob_threshold=1e-3, cooperativity=0.92):
    """
    Performs bounded multipath search, groups paths by terminal endpoint
    and length, and returns the partition function sum (Z) across degenerate paths.
    """
    import heapq
    from collections import defaultdict

    pq = []
    visited = set()

    # Store group data: endpoint -> list of (weight, path)
    endpoint_groups = defaultdict(list)

    for root in root_indices:
        if root in g:
            heapq.heappush(pq, (0.0, root, 0, [root]))

    while pq:
        curr_weight, u_node, depth, path = heapq.heappop(pq)

        state = (u_node, frozenset(path))
        if state in visited:
            continue
        visited.add(state)

        if len(path) > 1:
            prob = np.exp(-curr_weight)
            if prob >= prob_threshold:
                endpoint_groups[(u_node, len(path))].append((curr_weight, path))
                endpoint_groups[u_node].append((curr_weight, path))

        if depth >= max_depth:
            continue

        for v_node in g.neighbors(u_node):
            if v_node in path:
                continue

            edge_weight = g[u_node][v_node]['weight']
            next_weight = curr_weight + edge_weight * (cooperativity ** depth)
            next_prob = np.exp(-next_weight)

            if next_prob >= prob_threshold:
                next_path = path + [v_node]
                heapq.heappush(pq, (next_weight, v_node, depth + 1, next_path))

    # Compile partition sums and representative paths
    final_results = []
    for endpoint, paths_data in endpoint_groups.items():
        # Calculate Partition Sum: Z = sum(exp(-W_i))
        z_total = sum(np.exp(-w) for w, p in paths_data)

        # Representative path (lowest weight / highest probability)
        best_path = min(paths_data, key=lambda x: x[0])[1]

        final_results.append((best_path, float(z_total)))

    return final_results
