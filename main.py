# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    # pass

    sub_matrix = "./substitution_matrices/BLOSUM62.mat"
    
    nw_align = NeedlemanWunsch(sub_matrix, gap_open=-10, gap_extend=-1)

    #we calculate the alignment score for each of the species, including human

    hs_score, _, _ = nw_align.align(hs_seq, hs_seq)
    g_score, _, _ = nw_align.align(hs_seq, gg_seq)
    m_score, _, _ = nw_align.align(hs_seq, mm_seq)
    br_score, _, _ = nw_align.align(br_seq, gg_seq)
    tt_score, _, _ = nw_align.align(tt_seq, hs_seq)

    #make a dictionary of the scores and their species

    myscores = {'Homo_sapiens': hs_score, 'Gallus_gallus': g_score, 'Mus_musculus': m_score, 'Balaeniceps_rex': br_score, 'tursiops_truncatus': tt_score}

    #we print in order (greatest to least) the scores and the species
    species_list = []
    scores_list = []
    for w in sorted(myscores, key=myscores.get, reverse=True):
        species_list.append(w)
        scores_list.append(myscores[w])

    
    print(f"species order: {species_list}")
    print(f"scores of the species order: {scores_list}")

if __name__ == "__main__":
    main()
