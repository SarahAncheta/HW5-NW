# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing
        
        #we create empty matrices for the Ix, Iy, M and backtracking and after filling them in we store them as attributes of self (at the end)
        j_s = len(seqA) + 1
        i_s= len(seqB) + 1
        Ix = np.full((i_s, j_s), -np.inf) #the three matrices correspond to the three possible moves
        Iy = np.full((i_s, j_s), -np.inf)
        M = np.full((i_s, j_s), -np.inf)

        h = self.gap_open 
        g = self.gap_extend

       
        #points to M are stored as 0, Ix are 1, Iy are 2. when it has no pointer, it is stored as a -1
        backmatrix = np.full((i_s, j_s), -1)


        #we initialize the first row and column for each of the matrices

        for i in range(i_s): #initialize first column for the three matrices
            Iy[i][0] = h + g*i  # set first column value for Iy based on initial starting condition
            Ix[i][0]= -np.inf # first column of Ix is negative inf
            if i == 0:  # Set 0,0 item of M to 0, otherwise column value is negative inf
                M[i][0] = 0
                
            else:
                M[i][0] = -np.inf
                backmatrix[i][0] = 2 #points back to the Iy value

        for j in range(j_s): # initialize first row for the three matrices
            Ix[0][j] = h + g*j # first row value of Ix based on intitial starting condition
            Iy[0][j] = -np.inf # otherwise for Iy is negative inf
            if j == 0:        #same for M (neg inf), except for the 0,0 item which is 0
                M[0][j] = 0
                
            else:
                M[0][j] = -np.inf
                backmatrix[0][j] = 1 #points back to Ix value

        # TODO: Implement global alignment here    


        for i, j in np.ndindex(i_s, j_s): #we only iterate through the inside, ignore edges
            if i == 0 or j == 0:
                continue

            s_value = self.sub_dict[(seqB[i-1], seqA[j-1])] #in the following lines, we calculate the values of M, Ix and Iy

            m_vals = np.array([
                M[i-1][j-1] + s_value,
                Ix[i-1][j-1] + s_value,
                Iy[i-1][j-1] + s_value
            ])

            M[i][j]= np.max(m_vals)

            x_vals = np.array([
                M[i-1][j] + h + g,
                Ix[i-1][j] + g,
                -np.inf
            ])

            Ix[i][j] = np.max(x_vals)

            y_vals = np.array([
                M[i][j-1] + h + g,
                -np.inf,
                Iy[i][j-1] + g
            ])

            Iy[i][j] = np.max(y_vals)

            #our pointer stores the "point" back to the largest score
            pointer = np.argmax([M[i][j], Ix[i][j], Iy[i][j]])
            
            #we say which matrix it came from (0 indicates M, 1 indices Ix, 2 indicates Iy)
            if pointer == 0:
                backmatrix[i][j] = 0
            elif pointer == 1:
                backmatrix[i][j] = 1
            elif pointer == 2:
                backmatrix[i][j] = 2
            else:
                print('trap')

        #we store all these matrices we calculated as attributes

        self.Ix = Ix
        self.Iy = Iy
        self.M = M
        self.backmatrix = backmatrix

        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        max_row, max_col = self.M.shape
        max_row -= 1 # we subtract 1 because we have 0 indexing (starts with 0)
        max_col -= 1

        seqA = self._seqA
        seqB = self._seqB

        final_scores = np.array([self.M[max_row][max_col], self.Ix[max_row][max_col], self.Iy[max_row][max_col]])
        self.alignment_score = np.max(final_scores) #get the highest score across all three matrices, for the last/corner value

        row_val = max_row # we iterate backwards through the rows and columns
        col_val = max_col

        seqA_align = ''
        seqB_align = ''

        while row_val >= 0 or col_val >= 0: #checking if we hit the top right corner

            my_backtrack_index = self.backmatrix[row_val][col_val] #we get the pointer from the backmatrix for that particular square

            if my_backtrack_index == 0: #this is an M movement, take the diagonal

                seqA_align = seqA[col_val - 1] + seqA_align
                seqB_align = seqB[row_val - 1 ] + seqB_align

                row_val -= 1
                col_val -= 1

            if my_backtrack_index == 1: #this is an Ix movement, make a space in seqA and traverse up a row

                seqA_align = '-' + seqA_align
                seqB_align = seqB[row_val - 1] + seqB_align
                row_val -= 1


            if my_backtrack_index == 2: #this is an Iy movement, make a space in seqB and traverse left a column

                seqA_align = seqA[col_val - 1] + seqA_align
                seqB_align = '-' + seqB_align
                col_val -= 1

            if my_backtrack_index == -1: #if we get a -1 value, that means we have hit the edge
                if col_val == 0 and row_val == 0: # until we hit the top left corner, add spaces as needed
                    break

                if row_val == 0:
                    seqA_align = seqA[col_val - 1] + seqA_align
                    seqB_align = '-' + seqB_align
                    col_val -= 1

                if col_val == 0:
                    seqA_align = '-' + seqA_align
                    seqB_align = seqB[row_val - 1] + seqB_align
                    row_val -= 1

            
        self.seqA_align = seqA_align
        self.seqB_align = seqB_align

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header
