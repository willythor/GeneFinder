# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Willem Thorbecke

"""


import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq
dna = load_seq("./data/X73525.fa")

def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """

    if nucleotide == 'T':
        return 'A'
    elif nucleotide == 'A':
        return 'T'
    elif nucleotide == 'G':
        return 'C'
    else: 
        return 'G'
 


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    newDNA = ''
    length = len(dna)
    reverseDNA = dna[::-1]
   
    for x in range(0, length):
        newDNA = newDNA + get_complement(reverseDNA[x])
    return newDNA


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    x = 3 #starts searching after start sequence
    length = len(dna)
    while x < length:
        if dna[x:(x+3)] == 'TAG' or dna[x:(x+3)] == 'TAA' or dna[x:(x+3)] == 'TGA': 
            return dna[:x]
        else:
            x += 3
    return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    length = len(dna)
    x = 0
    non_nested_ORF = []
    while x < length:
        if dna[x:x+3] == 'ATG':
            section = rest_of_ORF(dna[x:])
            non_nested_ORF.append(section)
            x += len(section)
        else:    
            x += 3
    return non_nested_ORF 


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    ORF_Stag = [] 
    ORF_Stag = ORF_Stag + find_all_ORFs_oneframe(dna) #searches for ORFs with no stagger
    ORF_Stag = ORF_Stag + find_all_ORFs_oneframe(dna[1:]) #searches for ORFs with a one letter stagger
    ORF_Stag = ORF_Stag + find_all_ORFs_oneframe(dna[2:]) #searches for ORFs with a two letter stagger
    return ORF_Stag

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """

    R_Complements = get_reverse_complement(dna)
    ORFs = find_all_ORFs(dna) + find_all_ORFs(R_Complements)
    return ORFs

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    all_Orfs = find_all_ORFs_both_strands(dna)
    longest = 0
    length = len(all_Orfs)
    for x in range(0, length):
        if longest < len(all_Orfs[x]):
            longest = len(all_Orfs[x])
            position = x
        else:
            continue
    return all_Orfs[position]


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 
    """

    longest = 0

    dna_shuffle = shuffle_string(dna)

    for i in range(num_trials):
        dna_shuffle = shuffle_string(dna)

        longest_ORF_length = len(longest_ORF(dna_shuffle))

        if longest < longest_ORF_length:
            longest = longest_ORF_length
        else:
            continue
    return longest


def coding_strand_to_AA(dna):


    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents a protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
        >>> coding_strand_to_AA("TTTATCATGTTAGTTA")
        'FIMLV'
    """

    amino_acids_sequence = []
    x = 0
    while x < (len(dna)):
        codon = str(dna[x:x+3])
        if len(codon) == 3: #converts the gene sequence only if its 3 characters long
            holder = aa_table[codon]
            amino_acids_sequence.append(holder)
        else:
            return ''.join(amino_acids_sequence)
        x += 3
    return ''.join(amino_acids_sequence)


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """

    aa_Sequences = []

    threshold = longest_ORF_noncoding(dna, 1500)

    all_Orfs = find_all_ORFs_both_strands(dna)
    for orf in all_Orfs:
        if len(orf) > threshold:
            aa_Sequences.append(coding_strand_to_AA(orf))
        else:
            continue
    return aa_Sequences


print gene_finder(dna)

"""
    if __name__ == "__main__":  
    import doctest
    doctest.testmod()
    doctest.run_docstring_examples(coding_strand_to_AA, globals())
"""