import math
import random
import protein

def is_valid_protein_sequence(seq: str) -> bool:
    """ Check if the given sequence is a valid protein sequence.
        A valid protein sequence contains only the 20 standard amino acids represented by their single-letter codes.

        Args:
            seq (str): The input amino acid sequence.
            Sequence containing the 20 different amino-acids (RNDQEHKSTACGILMFPWYV), they must be upper case letters.

        Returns:
            bool: True if the sequence is valid, False otherwise.
    """
    valid_amino_acids = set('RNDQEHKSTACGILMFPWYV')
    for aa in seq:
        if aa not in valid_amino_acids:
            return False
    return True

def hp_seq (seq):
    """ Convert an amino acid sequence into a hydrophobic-polar (HP) sequence.
        Hydrophobic residues are represented by 'H' and polar residues by 'P'.

        Args:
            seq (str): The input amino acid sequence.
            Sequence containing the 20 different amino-acids (RNDQEHKSTACGILMFPWYV), they must be upper case letters.

        Returns:
            str: The corresponding HP sequence.
    """
    if not is_valid_protein_sequence(seq):
        raise ValueError("Invalid protein sequence. Only standard amino acid single-letter codes are allowed.")

    # Define hydrophobic and polar amino acids
    hydrophobic = set('AILMFWYV') # Hydrophobic amino acids
    polar = set('RNDQEHKSTC') # Polar amino acids

    # Convert the sequence to HP representation
    hp_sequence = ''.join(['H' if aa in hydrophobic else 'P' if aa in polar else 'X' for aa in seq])

    return hp_sequence

def is_valid_hp_sequence(hp_seq: str) -> bool:
    """ Check if the given sequence is a valid hydrophobic-polar (HP) sequence.
        A valid HP sequence contains only the characters 'H' and 'P'.

        Args:
            hp_seq (str): The input HP sequence.

        Returns:
            bool: True if the sequence is valid, False otherwise.
    """
    valid_chars = set('HP')
    for char in hp_seq:
        if char not in valid_chars:
            return False
    return True

def linear_structure(seq:str):
    """Generate a linear structure representation for a given HP sequence.

       Args:
           seq (str): The input HP sequence.

       Returns:
           dict: A dictionary mapping position indices to tuples containing (amino_acid, coordinates).
                Format: {0: ('H', (0, 0, 0)), 1: ('P', (1, 0, 0)), ...}
    """
    if not is_valid_hp_sequence(seq):
        raise ValueError("Invalid HP sequence. Only 'H' and 'P' characters are allowed.")

    structure = {}

    # Generate linear coordinates along the x-axis
    for i, residue in enumerate(seq):
        # Map position index to (amino_acid, coordinates) tuple
        structure[i] = (residue, (i, 0, 0))

    return structure

def get_distance(coord1: list, coord2: list) -> float:
    """
    Compute Euclidean distance between two 3D coordinates.
    
    Args:
        coord1 (list): list of 3 integers representing the first coordinates (x1, y1, z1).
        coord2 (list): list of 3 integers representing the second coordinates (x2, y2, z2).
        
    Returns:
        float: Euclidean distance between the two coordinates.
    """
    
    # Extract coordinates
    x1, y1, z1 = coord1
    x2, y2, z2 = coord2
    
    # Calculate Euclidean distance using the formula: sqrt((x2-x1)² + (y2-y1)² + (z2-z1)²)
    distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    
    return distance

def get_neighbors(hp_structure: dict, position:int) -> list:
    """
    Get which are the neighbors of a given amino acid in the HP structure.
    Args:
        hp_structure (dict): A dictionary mapping position indices to tuples containing (amino_acid, coordinates).
                             Format: {0: ('H', (0, 0, 0)), 1: ('P', (1, 0, 0)), ...}
        position (int): The index of the amino acid whose neighbors are to be found.

    Returns:
        list: A list of indices representing the neighboring amino acids.
    """
    if position not in hp_structure:
        raise ValueError("Position not found in the HP structure.")
    
    neighbors = []
    target_coord = hp_structure[position][1]
    
    for aa, coord in hp_structure.items():
        if aa != position and aa != position - 1 and aa != position + 1:  # Exclude self and sequential neighbors
            distance = get_distance(target_coord, coord[1])
            if abs(distance - 1.0) < 1e-10:  # Only exactly 1 unit apart (not diagonal)
                neighbors.append(hp_structure[aa][0])

    return neighbors

def energy(hp_structure: dict, e=1.0) -> float:
    """
    Calculate the energy of a given HP structure based on non-sequential H-H contacts.
    
    Args:
        hp_structure (dict): A dictionary mapping position indices to tuples containing (amino_acid, coordinates).
                             Format: {0: ('H', (0, 0, 0)), 1: ('P', (1, 0, 0)), ...}
        e (float): The energy contribution for each non-sequential H-H contact. Default is 1.0.
                                                  
    Returns:
        float: The calculated energy of the structure.
    """
    energy = 0
    h_h_counts = 0
    num_aa = len(hp_structure)
    
    for i in range(num_aa):
        aa_i, coord_i = hp_structure[i]
        if aa_i == 'H':
            for j in get_neighbors(hp_structure, i):
                if j == 'H':
                    h_h_counts += 1
    
    energy = - (h_h_counts / 2) * e  # Each contact counted twice

    return energy