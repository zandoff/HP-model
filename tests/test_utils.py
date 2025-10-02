import unittest
import math
from utils import (
    is_valid_protein_sequence, 
    hp_seq, 
    is_valid_hp_sequence, 
    linear_structure, 
    get_distance, 
    get_neighbors, 
    energy, 
    folding
)


class TestUtils(unittest.TestCase):

    def test_is_valid_protein_sequence(self):
        """Test protein sequence validation."""
        # Valid sequences
        self.assertTrue(is_valid_protein_sequence("ACDEFGHIKLMNPQRSTVWY"))
        self.assertTrue(is_valid_protein_sequence("ALANINE"))
        self.assertTrue(is_valid_protein_sequence(""))
        
        # Invalid sequences
        self.assertFalse(is_valid_protein_sequence("ACDEFGHIKLMNPQRSTVWYX"))  # Contains X
        self.assertFalse(is_valid_protein_sequence("acdef"))  # Lowercase
        self.assertFalse(is_valid_protein_sequence("ALAN1NE"))  # Contains number
        self.assertFalse(is_valid_protein_sequence("ALAN-INE"))  # Contains hyphen

    def test_hp_seq(self):
        """Test amino acid to HP sequence conversion."""
        # Test hydrophobic amino acids
        self.assertEqual(hp_seq("AILMFWYV"), "HHHHHHHH")
        
        # Test polar amino acids  
        self.assertEqual(hp_seq("RNDQEHKSTC"), "PPPPPPPPPP")
        
        # Test mixed sequence
        self.assertEqual(hp_seq("ALANINE"), "HPHPHHP")
        
        # Test empty sequence
        self.assertEqual(hp_seq(""), "")
        
        # Test invalid sequence should raise error
        with self.assertRaises(ValueError):
            hp_seq("ALANX")

    def test_is_valid_hp_sequence(self):
        """Test HP sequence validation."""
        # Valid HP sequences
        self.assertTrue(is_valid_hp_sequence("HPHP"))
        self.assertTrue(is_valid_hp_sequence("HHHHPPPP"))
        self.assertTrue(is_valid_hp_sequence(""))
        self.assertTrue(is_valid_hp_sequence("H"))
        self.assertTrue(is_valid_hp_sequence("P"))
        
        # Invalid HP sequences
        self.assertFalse(is_valid_hp_sequence("HPXA"))
        self.assertFalse(is_valid_hp_sequence("hphp"))  # Lowercase
        self.assertFalse(is_valid_hp_sequence("HP1P"))  # Contains number

    def test_linear_structure(self):
        """Test linear structure generation."""
        # Test simple sequence
        result = linear_structure("HP")
        expected = {0: ('H', (0, 0, 0)), 1: ('P', (1, 0, 0))}
        self.assertEqual(result, expected)
        
        # Test longer sequence
        result = linear_structure("HPHP")
        expected = {
            0: ('H', (0, 0, 0)), 
            1: ('P', (1, 0, 0)), 
            2: ('H', (2, 0, 0)), 
            3: ('P', (3, 0, 0))
        }
        self.assertEqual(result, expected)
        
        # Test empty sequence
        result = linear_structure("")
        self.assertEqual(result, {})
        
        # Test invalid sequence should raise error
        with self.assertRaises(ValueError):
            linear_structure("HPX")

    def test_get_distance(self):
        """Test Euclidean distance calculation."""
        # Test same point
        self.assertEqual(get_distance([0, 0, 0], [0, 0, 0]), 0.0)
        
        # Test unit distances
        self.assertEqual(get_distance([0, 0, 0], [1, 0, 0]), 1.0)
        self.assertEqual(get_distance([0, 0, 0], [0, 1, 0]), 1.0)
        self.assertEqual(get_distance([0, 0, 0], [0, 0, 1]), 1.0)
        
        # Test diagonal distance
        self.assertAlmostEqual(get_distance([0, 0, 0], [1, 1, 0]), math.sqrt(2), places=10)
        
        # Test 3-4-5 triangle
        self.assertEqual(get_distance([0, 0, 0], [3, 4, 0]), 5.0)
        
        # Test negative coordinates
        self.assertEqual(get_distance([-1, -1, -1], [1, 1, 1]), math.sqrt(12))

    def test_get_neighbors(self):
        """Test neighbor detection."""
        # Test structure with no spatial neighbors
        structure = {
            0: ('H', (0, 0, 0)),
            1: ('P', (1, 0, 0)),
            2: ('H', (2, 0, 0))
        }
        self.assertEqual(get_neighbors(structure, 1), [])
        
        # Test structure with spatial neighbors
        structure = {
            0: ('H', (0, 0, 0)),
            1: ('P', (1, 0, 0)),    # Sequential neighbor - excluded
            2: ('H', (0, 1, 0)),    # Spatial neighbor at distance 1.0
            3: ('P', (1, 1, 0)),    # Not a neighbor of position 0 (distance sqrt(2))
            4: ('H', (0, 0, 1))     # Spatial neighbor at distance 1.0
        }
        neighbors = get_neighbors(structure, 0)
        self.assertIn('H', neighbors)  # From position 2
        self.assertIn('H', neighbors)  # From position 4
        self.assertEqual(len(neighbors), 2)
        
        # Test invalid position
        with self.assertRaises(ValueError):
            get_neighbors(structure, 10)

    def test_energy(self):
        """Test energy calculation."""
        # Test structure with no H-H contacts
        structure = {
            0: ('H', (0, 0, 0)),
            1: ('P', (1, 0, 0)),
            2: ('H', (2, 0, 0))
        }
        self.assertEqual(energy(structure), 0.0)
        
        # Test structure with H-H contacts
        structure = {
            0: ('H', (0, 0, 0)),
            1: ('P', (1, 0, 0)),    # Sequential neighbor - excluded
            2: ('H', (0, 1, 0)),    # H-H contact with position 0
            3: ('P', (2, 0, 0))
        }
        # Should have 1 H-H contact, so energy = -1.0
        self.assertEqual(energy(structure), -1.0)
        
        # Test with custom energy parameter
        self.assertEqual(energy(structure, e=2.0), -2.0)
        
        # Test structure with multiple H-H contacts
        structure = {
            0: ('H', (0, 0, 0)),
            1: ('P', (1, 0, 0)),
            2: ('H', (0, 1, 0)),    # H-H contact with position 0
            3: ('H', (0, 0, 1)),    # H-H contact with position 0
            4: ('P', (2, 0, 0))
        }
        # Should have 2 H-H contacts, so energy = -2.0
        self.assertEqual(energy(structure), -2.0)

    def test_folding(self):
        """Test folding function (placeholder)."""
        # Test that folding returns the same structure (placeholder behavior)
        structure = {
            0: ('H', (0, 0, 0)),
            1: ('P', (1, 0, 0)),
            2: ('H', (2, 0, 0))
        }
        result = folding(structure)
        self.assertEqual(result, structure)
        
        # Test with different max_iterations parameter
        result = folding(structure, max_iterations=500)
        self.assertEqual(result, structure)

    def test_integration(self):
        """Test integration of multiple functions."""
        # Create a complete workflow
        protein_seq = "ALAMINA"
        hp_sequence = hp_seq(protein_seq)
        self.assertTrue(is_valid_hp_sequence(hp_sequence))
        
        structure = linear_structure(hp_sequence)
        self.assertEqual(len(structure), len(hp_sequence))
        
        # Calculate energy of linear structure
        initial_energy = energy(structure)
        
        # Fold structure (placeholder)
        folded_structure = folding(structure)
        final_energy = energy(folded_structure)
        
        # Energy should be the same since folding is a placeholder
        self.assertEqual(initial_energy, final_energy)


if __name__ == '__main__':
    # Run the tests
    unittest.main(verbosity=2)