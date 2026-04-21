import unittest
import os
import tempfile
import networkx as nx
import MDAnalysis as mda
import numpy as np

from water_bridges_nw.math_utils import switching_function, calculate_hbond_probability
from water_bridges_nw.core import build_graph, compute_edge_probabilities, traverse_network
from unittest.mock import patch
from water_bridges_nw.visualize import export_vmd_script, export_pymol_script, run_visualization

class TestMathUtils(unittest.TestCase):
    def test_switching_function(self):
        # When dist == threshold, it should return the limit (power_num / power_den)
        # For powers 8 and 12, limit is 8/12 = 2/3
        self.assertAlmostEqual(switching_function(1.0, 1.0, 8, 12), 2.0 / 3.0)

        # When dist < threshold, it should be high
        prob_high = switching_function(0.5, 1.0)
        self.assertTrue(0 < prob_high <= 1)

    def test_calculate_hbond_probability(self):
        # Ideal mock h-bond values
        prob = calculate_hbond_probability(2.7, 1.0, 1.7)
        self.assertTrue(prob >= 0.0)

class TestCoreFunctions(unittest.TestCase):
    def setUp(self):
        # Create a tiny mock topology in memory using MDAnalysis Universe
        # 3 waters, 1 ligand
        self.u = mda.Universe.empty(4, trajectory=True, n_residues=4, atom_resindex=[0, 1, 2, 3])
        self.u.add_TopologyAttr('resname', ['LIG', 'SOL', 'SOL', 'SOL'])
        self.u.add_TopologyAttr('name', ['O1', 'OW', 'OW', 'OW'])
        self.u.atoms.positions = np.array([
            [0.0, 0.0, 0.0],
            [2.5, 0.0, 0.0], # Connected to LIG
            [5.0, 0.0, 0.0], # Connected to first water
            [10.0, 0.0, 0.0] # Far away, no connection
        ])
        # Box dimensions needed for capped_distance
        self.u.dimensions = np.array([100.0, 100.0, 100.0, 90.0, 90.0, 90.0])

    def test_build_graph(self):
        water_atoms = self.u.select_atoms("resname SOL")
        root_atoms = self.u.select_atoms("resname LIG")
        g, roots = build_graph(self.u, water_atoms, root_atoms, max_distance=3.5)
        self.assertEqual(len(roots), 1)
        self.assertEqual(roots[0], 0)
        # Edges expected: 0-1 (dist 2.5), 1-2 (dist 2.5)
        self.assertIn((0, 1), g.edges)
        self.assertIn((1, 2), g.edges)
        self.assertNotIn((2, 3), g.edges) # Distance 5.0 is > 3.5

    def test_traverse_network(self):
        g = nx.Graph()
        g.add_node(0)
        g.add_node(1)
        g.add_node(2)

        g.add_edge(0, 1, prob=0.8, weight=-np.log(0.8))
        g.add_edge(1, 2, prob=0.8, weight=-np.log(0.8))

        paths = traverse_network(g, [0], max_depth=5, prob_threshold=0.5)

        # Valid paths should include [0, 1] and [0, 1, 2]
        # (Actually, because of the loop structure, when it expands to 2,
        # [0, 1] might just expand and not be recorded separately if we only record terminal paths,
        # but let's just check length)

        self.assertTrue(len(paths) > 0)

        # We should find the path to node 2
        found_longest = False
        for p, prob in paths:
            if p == [0, 1, 2]:
                found_longest = True
                self.assertAlmostEqual(prob, 0.64)
        self.assertTrue(found_longest)

class TestVisualization(unittest.TestCase):
    def test_export_vmd_script(self):
        # We need to mock a jsonl file now since visualization reads files directly
        with tempfile.NamedTemporaryFile(delete=False, suffix=".jsonl", mode='w') as f_json:
            json_name = f_json.name
            # Write mock frame data
            f_json.write('{"type": "metadata"}\n')
            f_json.write('{"type": "frame", "frame_idx": 0, "paths": [{"nodes": [0, 1], "coords": [[0,0,0], [1,1,1]]}]}\n')

        with tempfile.NamedTemporaryFile(delete=False, suffix=".tcl") as f:
            temp_name = f.name

        export_vmd_script(json_name, output_file=temp_name, mode="frame", frame_idx=0)

        with open(temp_name, 'r') as f:
            content = f.read()
            self.assertIn("mol selection \"index 0 1\"", content)

        os.remove(temp_name)
        os.remove(json_name)

    @patch('water_bridges_nw.visualize.logger.error')
    def test_run_visualization_file_not_found(self, mock_logger_error):
        run_visualization("non_existent_file.jsonl", format="vmd")
        mock_logger_error.assert_called_once_with("Data file non_existent_file.jsonl not found.")

    @patch('water_bridges_nw.visualize.export_vmd_script')
    def test_run_visualization_vmd_routing(self, mock_export):
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            pass
        try:
            run_visualization(temp_file.name, format="vmd", mode="frame", frame_idx=0, output_file="test_out")
            mock_export.assert_called_once_with(temp_file.name, output_file="test_out.tcl", mode="frame", frame_idx=0)
        finally:
            os.remove(temp_file.name)

    @patch('water_bridges_nw.visualize.export_pymol_script')
    def test_run_visualization_pymol_routing(self, mock_export):
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            pass
        try:
            run_visualization(temp_file.name, format="pymol", mode="density", output_file="test_out")
            mock_export.assert_called_once_with(temp_file.name, output_file="test_out.py", mode="density", frame_idx=None)
        finally:
            os.remove(temp_file.name)

    @patch('water_bridges_nw.visualize.export_chimera_script')
    def test_run_visualization_chimera_routing(self, mock_export):
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            pass
        try:
            run_visualization(temp_file.name, format="chimera", mode="frame", frame_idx=10, output_file="test_out.py")
            mock_export.assert_called_once_with(temp_file.name, output_file="test_out.py", mode="frame", frame_idx=10)
        finally:
            os.remove(temp_file.name)

    @patch('water_bridges_nw.visualize.logger.error')
    def test_run_visualization_unknown_format(self, mock_logger_error):
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            pass
        try:
            run_visualization(temp_file.name, format="unknown", mode="density", output_file="test_out")
            mock_logger_error.assert_called_once_with("Unknown format: unknown")
        finally:
            os.remove(temp_file.name)

if __name__ == '__main__':
    unittest.main()
