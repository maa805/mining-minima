class MiningMinima:
    '''
    A single MM style simulation
    '''

    def __init__(self, seq1='', seq2='', infile='')

    '''
    Create a MM object.
    
    Parameters
    ---------------------------
    seq1:   The sequence for a single strand simulation or the first strand in a double strand simulation
    seq2:   The second sequence for a double strand simulation
    infile: A PDB file on which to run mining minima
    '''
    if !infile: 

        self.seq1 = seq1
        if seq2: self.seq2 = seq2

    else: self.infile = infile
    
