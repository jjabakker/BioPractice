from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.Data import CodonTable
from Bio import SeqIO
from Bio import Align
import Bio
import pprint


def pattern_count(seq_to_search, pattern):
    # myseq: Seq = Seq (seq_to_search, IUPAC.unambiguous_dna)

    positions = []
    pos = seq_to_search.find(pattern)
    while pos != -1:
        positions.append(pos)
        pos = seq_to_search.find(pattern, pos + 1)
    return len(positions), positions


def frequent_words(seq_to_search, kmer):
    """

    Find the most frequent occurring k-mer, using a naive approach

    :param seq_to_search: the string in which you are searching
    :param kmer: the length of the repeat element
    :return:
        maxcount: the number of times the kmer was encountered
        str:      the kmer that was found
        position: the

    """

    startindex = 0
    maxcount = 0
    strings_searched = set()

    result_index = [0] * len(seq_to_search)
    for startindex in range(0, len(seq_to_search) - kmer + 1):
        str_search = seq_to_search[startindex:startindex + kmer]

        # If you have already searched for this string, do not do it again

        if str_search in strings_searched:
            continue
        else:
            strings_searched.add(str_search)
            nr_found, positions = pattern_count(seq_to_search, seq_to_search[startindex:startindex + kmer])
            if nr_found > maxcount:
                max_strings_found = []
                max_strings_found.append([seq_to_search[startindex:startindex + kmer], positions])

                maxcount = nr_found
                index = startindex
            elif nr_found == maxcount:
                max_strings_found.append([seq_to_search[startindex:startindex + kmer], positions])
            for i in positions:
                result_index[i] = nr_found

    return maxcount, max_strings_found


def wrapper_for_find_most_occurring_kmer(sequence, kmer):
    (n, seqs_found) = frequent_words(sequence, kmer)
    print('\n')
    if n > 1:
        for seq_found in seqs_found:
            print(f"String {seq_found[0]} with length {len(seq_found[0])} occurs {n} times op positions {seq_found[1]}")
            i = seq_found[1][0]
            i = seq_found[1][1]
            print_sequence(sequence, 100)
            print('\n')

            for i in range(0, n):
                str = sequence[seq_found[1][i]:seq_found[1][i] + len(seq_found[0])]
                print(str)
            print('\n')

    else:
        print(f"There is no {kmer}-kmer")

    print('\n')


def find_highest_kmer_and_segment(seq, segment_len, kmer_length):

    """
    The sequence is cut in parts of seqment_length and in each part the highest repeat occurrence is determined
    :param segment_len:
    :param kmer_length:
    :return:
    """

    seqlen = len(seq)
    i = 0
    maxcount = 0
    maxindex = 0

    while i + segment_len < seqlen:
        count, max_strings_found = frequent_words(seq[i:i + segment_len], kmer_length)
        if count > maxcount:
            maxcount = count
            maxindex = i
            print(maxcount, maxindex)
        i += segment_len

    wrapper_for_find_most_occurring_kmer(seq[maxindex:maxindex + segment_len], kmer_length)


def find_highest_GC():

    """
    Find the AT percentages in segments of 500
    Print what the percentage is and where it wss found
    """

    ori = SeqIO.read('ori.fasta.txt', 'fasta')
    genome = SeqIO.read('vibrio_cholerae.fasta.txt', 'fasta')
    print(ori.seq)

    print(f'The AT percentage in ori is {100 - GC(ori.seq):5.2f}')
    print(f'The AT percentage in genome is {100 - GC(genome.seq):5.2f}')

    i = 0
    seq = genome.seq
    seqlen = len(genome.seq)
    maxperc = 0
    maxindex = 0

    while i + 500 < seqlen:
        perc = GC(seq[i:i + 500])
        if perc > maxperc:
            maxperc = perc
            maxindex = i
        i += 500
    print(maxperc, maxindex)


def print_sequence(seq, line_length):
    """
    Print a specified sequence in lines of line_lengtg
    :param seq:
    :param line_length:
    """

    i = 0
    header = ""
    while i < line_length / 10:
        header = header + 10*str(i)
        i += 1
    print(header)
    header = int(line_length / 10) * '0123456789'
    print(header+'\n')
    while len(seq) > line_length:
        print(seq[:line_length])
        seq = seq[i + line_length:]
    if len(seq) > 0:
        print(seq)


def analyse_ori(ori):

    wrapper_for_find_most_occurring_kmer(ori, 9)

    genome = SeqIO.read('vibrio_cholerae.fasta.txt', 'fasta')
    count, positions = pattern_count(genome.seq, "ATGATCAAG")
    count, positions = pattern_count(genome.seq, "CTTGATCAT")

    m = 1


def main():
    print(pattern_count("GACCATCAAAACTGATAAACTACTTAAAAATCAGT", "AAA"))
    print(pattern_count("GACCATCAAAACTGATAAACTACTTAAAAATCAGT", "GATA"))
    print(pattern_count("GACCATCAAAACTGATAAACTACTTAAAAATCAGT", "AA"))
    print('\n\n')

    #
    # Find the occurence of kmers is
    #

    ori = SeqIO.read('ori.fasta.txt', 'fasta')

    k = 3
    count, max_strings_found = frequent_words(str(ori.seq), k)
    print(f'The maximum time a {k}-kmer was found is {count}')
    print(f'The number of patterns for that maximum is {len(max_strings_found)}')
    print('The positions are:\n')
    for l in range(len(max_strings_found)):
        print(max_strings_found[l])
    print('\n')

    k = 5
    count, max_strings_found = frequent_words(str(ori.seq), k)
    print(f'The maximum time a {k}-kmer was found is {count}')
    print(f'The number of patterns for that maximum is {len(max_strings_found)}')
    print('The positions are:\n\n')
    for l in range(len(max_strings_found)):
        print(max_strings_found[l])
    print('\n')

    k = 7
    count, max_strings_found = frequent_words(str(ori.seq), k)
    print(f'The maximum time a {k}-kmer was found is {count}')
    print(f'The number of patterns for that maximum is {len(max_strings_found)}')
    print('The positions are:\n\n')
    for l in range(len(max_strings_found)):
        print(max_strings_found[l])
    print('\n')


    k = 9
    count, max_strings_found = frequent_words(str(ori.seq), k)
    print(f'The maximum time a {k}-kmer was found is {count}')
    print(f'The number of patterns for that maximum is {len(max_strings_found)}')
    print('The positions are:\n\n')
    for l in range(len(max_strings_found)):
        print(max_strings_found[l])
    print('\n')

    '''
    wrapper_for_find_most_occurring_kmer('TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT', 2)
    wrapper_for_find_most_occurring_kmer('TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT', 3)
    wrapper_for_find_most_occurring_kmer('TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT', 4)
    wrapper_for_find_most_occurring_kmer('TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT', 5)
    wrapper_for_find_most_occurring_kmer('TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT', 6)
    wrapper_for_find_most_occurring_kmer('TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT', 7)
    wrapper_for_find_most_occurring_kmer('TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT', 8)
    '''

    count, positions = pattern_count("GACGATATACGACGATA", "ATA")
    print(positions)

    myseq = Seq("CCAGATC", IUPAC.unambiguous_dna)
    print(f'The reverse complement of {myseq} is: {myseq.reverse_complement()}')

    find_highest_GC()

    genome = SeqIO.read('vibrio_cholerae.fasta.txt', 'fasta')
    seq = genome.seq

    print(f'Length of genome is {len(seq)}')

    # find_highest_kmer_and_segment(seq, 500, 5)
    # find_highest_kmer_and_segment(seq, 500, 8)
    # find_highest_kmer_and_segment(seq, 500, 12)
    # find_highest_kmer_and_segment(seq, 500, 14)
    # find_highest_kmer_and_segment(seq, 500, 16)
    # find_highest_kmer_and_segment(seq, 500, 18)
    # find_highest_kmer_and_segment(seq, 500, 20)
    # find_highest_kmer_and_segment(seq, 500, 24)
    # find_highest_kmer_and_segment(seq, 500, 26)
    # find_highest_kmer_and_segment(seq, 500, 28)

    ori = SeqIO.read('ori.fasta.txt', 'fasta')

    # print_sequence(ori.seq, 120)
    # print('\n\n')
    # print_sequence(seq[107300:107800], 100)

    analyse_ori(ori.seq)


if __name__ == "__main__":
    main()
