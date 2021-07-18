# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 12:57:21 2021

@author: ariyaths
"""

from Bio import SeqIO


def create_dna( num, alphabet='actg'):
    """
    num : integer
    alphabet : string, optional, The default is 'actg'.
    Returns : string

    To create a random string of dna of desired length.
    """
    import random

    return ''.join([random.choice(alphabet) for i in range(num)])


def explore_blast(dna):
    """
    dna : string
    Returns: None

    To find a dna pattern in NCBI database using Blast
    """
    from Bio.Blast import NCBIWWW, NCBIXML

    result_handle = NCBIWWW.qblast("blastn", "nt", dna)
    blast_record = NCBIXML.read(result_handle)
    print(dna)
    print(len(blast_record.alignments))
    e_value_thresh = 0.01

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < e_value_thresh:
                print('****Alignment****')
                print('sequence: ', alignment.title)
                print('length: ', alignment.length)
                print('e value: ', hsp.expect)
                print(hsp.query)
                print(hsp.match)


def print_fasta_info(rec_dict):
    """
    rec_dict : dict of multiple dna sequence names and sequences
    Returns: None

    To describe the contents of the multi-dna-dict
    """

    def maplen(dictrow):
        return(dictrow[0], len(dictrow[1]))

    print("\n TOTAL NUMBER OF RECORDS : ", len(rec_dict), "\n")

    gene_len = dict(map(maplen, rec_dict.items()))
    print(*gene_len.items(),'\n',sep='\n')

    max_gene_len = max(gene_len.values())
    min_gene_len = min(gene_len.values())
    print("MAX GENE LENGTH : ", max_gene_len)
    print("MIN GENE LENGTH : ", min_gene_len)

    max_count = 0
    min_count = 0

    for k, val in gene_len.items():
        if val == max_gene_len:
            max_count += 1
            print("MAX KEY : ", k)
        if val == min_gene_len:
            min_count += 1
            print("MIN KEY : ", k)

    print("MAX COUNT : ", max_count)
    print("MIN COUNT : ", min_count)
    print('\n')


def orf_lengths(r_dict, orf_start=0):
    """    These are called reading frames 1, 2, and 3 respectively.
    An open reading frame (ORF) is the part of a reading frame that has the
    potential to encode a protein. It starts with a start codon (ATG), and ends
    with a stop codon (TAA, TAG or TGA)

    USE THIS FOR TESTING
    dna="atgaatgaatgcctaatagtgaatgccctgaataatagtgaatgcccccctgaataatagtgaatgccccccccctga"
    """

    def next_start_codon(dna):

        for i in range(0,len(dna),3):
            codon = dna[i:i+3].lower()
            if codon == 'atg':
                return i
        return i

    def next_stop_codon(dna):

        stop_codons = ['tga', 'tag', 'taa']
        for i in range(0,len(dna),3):
            codon = dna[i:i+3].lower()
            if codon in stop_codons:
                return i+3
        return 0

    result_list = []

    for seq_id in r_dict:
        inner_dict = {} #contains orf as keys and length as values
        outer_dict = {} #contains seq_id as keys and inner_dict{} as values

        frame = orf_start # start of open reading frame
        pos = 0 + frame
        pos_start = 0
        orf_len = 0

        dna = str(r_dict[seq_id].seq)

        while pos < len(dna):
            pos_start = next_start_codon(dna[pos:])
            orf_len = next_stop_codon(dna[pos + pos_start:])

            if orf_len == 0:
                break

            inner_dict[dna[pos + pos_start : pos + pos_start + orf_len]] = orf_len
            pos = pos + pos_start + orf_len    # new start position

        outer_dict[seq_id] = inner_dict

        if inner_dict == {}:
            max_key =''; max_value=0; freq = 0
        else:
            max_key = max(inner_dict, key=inner_dict.get)
            max_value = max(inner_dict.values())
            freq = sum(x == max_value for x in inner_dict.values())

        result_list.append((max_value, seq_id, freq, \
                            dna.find(max_key)+1, max_key[:40]))

    return result_list


def loop_repeats(r_dict, num):
    """ A repeat is a substring of a DNA sequence that occurs in multiple
    copies (more than one) somewhere in the sequence. Although repeats can
    occur on both the forward and reverse strands of the DNA sequence, we
    will only consider repeats on the forward strand here. Also we will
    allow repeats to overlap themselves. For example, the sequence ACACA
    contains two copies of the sequence ACA - once at position 1
    (index 0 in Python), and once at position 3.
    Given a length n, your program should be able to identify all repeats
    of length n in all sequences in the FASTA file. Your program should
    also determine how many times each repeat occurs in the file, and which
    is the most frequent repeat of a given length.

    return repeat_list, repeat_frequency, max_repeat_seq

    Parameters
    ----------
    dnadata : str
        A dna sequence.
    k : int
        length if k-mer.
    sorted_list : bool, optional
        If sorted_list, output is sorted. The default is False.
    Returns
    -------
    dict
        dict of {k-mer:count}.
    """

    def repeats(dnadata, k_len):

        freqdict = {}

        for i in range(len(dnadata) -k_len +1):
            kmer = dnadata[i:i+k_len]

            try:
                freqdict[kmer] += 1
            except KeyError:
                freqdict[kmer] = 1

        # Removing all kmer repeats less than 1
        freqdict_r = {kmer:count for kmer,count in freqdict.items() if count > 1}

        return freqdict_r


    seqdict = {}

    for seq_id in r_dict:
        dna = str(r_dict[seq_id].seq)
        seqdict[seq_id] = repeats(dna, num)

    repeatdict = {}

    for seq_id in seqdict:
        if seqdict[seq_id]:
            for k, value in seqdict[seq_id].items():
                if k in repeatdict:
                    repeatdict[k] += value
                else:
                    repeatdict[k] = value

    return seqdict, repeatdict


def print_repeat_info(r_dict, num):
    """
    r_dict : dict of multiple dna sequence names and sequences
    Returns: None

    To analyze the repeats in a multi-dna-dict
    """

    seqdict, repeatdict = loop_repeats(r_dict, num)

    for seq_id in seqdict:
        # print(seq_id, *seqdict[seq_id].items(), sep='\n')
        leng = [k for k,v in seqdict[seq_id].items() \
            if v == max(seqdict[seq_id].values())]
        print ('MAX = ', max(seqdict[seq_id].values()) \
               if leng else 'None', ' ', leng,'\n')


    # for k,v in repeatdict.items():
    #     if v > 2: print (k, v)

    lis = [ke for ke,va in repeatdict.items() if va == max(repeatdict.values())]
    print ('\n', 'MAX =', max(repeatdict.values()), lis)


if __name__ == '__main__':

    # dna = create_dna(10000)
    # explore_blast(dna)

    # FASTA_FILE = "data//fasta//dna_example.fasta"
    FASTA_FILE = "data//fasta//dna2.fasta"
    record_dict = SeqIO.index(FASTA_FILE, "fasta")

    print_fasta_info(record_dict)
    print(*orf_lengths(record_dict, orf_start = 2),'\n',sep='\n')
    print_repeat_info(record_dict, 7)

    record_dict.close()
