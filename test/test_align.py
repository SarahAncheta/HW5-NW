# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np


def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    sub_matrix = "./substitution_matrices/BLOSUM62.mat"
    
    nw_align = NeedlemanWunsch(sub_matrix, gap_open=-10, gap_extend=-1)
    score, _, _ = nw_align.align(seq1, seq2)

    assert score == 4
    test_M_array = np.array([[0.0, -np.inf, -np.inf, -np.inf, -np.inf], 
                                   [-np.inf,   5, -12, -12, -14], 
                                   [-np.inf, -11,   4, -1, -6], 
                                   [-np.inf, -13, -8,   5,  4]])

    assert np.allclose(nw_align.M, test_M_array)

    test_Ix_array = np.array([[-10, -11, -12, -13, -14], 
               [-np.inf, -12, -13, -14, -15],
               [-np.inf,  -6, -14, -15, -16],
               [-np.inf,  -7,  -7, -12, -17]])
    
    assert np.allclose(nw_align.Ix, test_Ix_array)

    test_Iy_array = np.array([[-np.inf, -np.inf, -np.inf, -np.inf, -np.inf],
                              [-11, -12,  -6,  -7,  -8],
                              [-12, -13, -14,  -7,  -8],
                              [-13, -14, -15, -16,  -6]])
    
    assert np.allclose(nw_align.Iy, test_Iy_array)

    

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    sub_matrix = "./substitution_matrices/BLOSUM62.mat"
    
    nw_align = NeedlemanWunsch(sub_matrix, gap_open=-10, gap_extend=-1)
    score, align3, align4 = nw_align.align(seq3, seq4)

    assert align3 == 'MAVHQLIRRP'
    assert align4 == 'M---QLIRHP'
    assert score == 17




