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
    score, align1, align2 = nw_align.align(seq1, seq2)

    print(nw_align.Ix)
    # print(nw_align.Iy)
    # print(nw_align.M)
    print(seq1)
    print(seq2)
    print(align1)
    print(align2)

    

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
    pass




