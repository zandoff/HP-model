def hp_seq (seq):
    """ Convert an amino acid sequence into a hydrophobic-polar (HP) sequence.
        Hydrophobic residues are represented by 'H' and polar residues by 'P'.
        
        Args:
            seq (str): The input amino acid sequence.
            Sequence containing the 20 different amino-acids (RNDQEHKSTACGILMFPWYV), they must be upper case letters.
        
        Returns:
            str: The corresponding HP sequence.
    """
    # Define hydrophobic and polar amino acids
    hydrophobic = set('AILMFWYV') # Hydrophobic amino acids
    polar = set('RNDQEHKSTC') # Polar amino acids

    # Convert the sequence to HP representation
    hp_sequence = ''.join(['H' if aa in hydrophobic else 'P' if aa in polar else 'X' for aa in seq])

    return hp_sequence

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
           list: A list of tuples representing the coordinates of each residue in a linear structure.
    """
       