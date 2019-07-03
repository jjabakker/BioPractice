from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.Data import CodonTable
from Bio import SeqIO
from Bio import Align
import Bio
import pprint

def print_header(page, text=''):

    print('\n')
    print('-' * 80)
    if page != 0:
        print(f'\tProgramming - Page: {page}  {text}')
    elif text != '':
        print(f'\tProgramming: {text}')
    else:
        print('\n')
    print('-' * 80)
    print('\n')


def pattern_count(seq_to_search, pattern):
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

    max_count = 0
    strings_searched = set()
    max_strings_found = []

    result_index = [0] * len(seq_to_search)
    for startindex in range(0, len(seq_to_search) - kmer + 1):
        str_search = seq_to_search[startindex:startindex + kmer]

        # If you have already searched for this string, do not do it again
        if str_search in strings_searched:
            continue
        else:
            strings_searched.add(str_search)
            nr_found, positions = pattern_count(seq_to_search, seq_to_search[startindex:startindex + kmer])
            if nr_found > 1:
                if nr_found > max_count:
                    max_strings_found.append([seq_to_search[startindex:startindex + kmer], positions])
                    max_count = nr_found
                    index = startindex
                elif nr_found == max_count:
                    max_strings_found.append([seq_to_search[startindex:startindex + kmer], positions])
                for i in positions:
                    result_index[i] = nr_found

    return max_count, max_strings_found


'''
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

    print_frequent_words(seq[maxindex:maxindex + segment_len], kmer_length)
'''


def find_max_GC(sequence: object, block_size: object = 0):
    """
    Find the AT percentages in segments of 500
    Print what the percentage is and where it wss found
    """

    if block_size == 0:
        block_size = len(sequence)

    i = 0
    seqlen = len(sequence)
    maxperc = 0
    maxindex = 0

    while i + block_size <= seqlen:
        perc = GC(sequence[i:i + block_size])
        if perc > maxperc:
            maxperc = perc
            maxindex = i
        i += block_size
    return maxperc, maxindex


def print_sequence(seq, line_length, has_header = False):

    i = 0
    if has_header:
        header = ""
        while i < line_length / 10:
            header = header + 10 * str(i)
            i += 1
        print(header)
        header = int(line_length / 10) * '0123456789'
        print(header)
        print('-'*line_length)

    while len(seq) > line_length:
        print(seq[:line_length])
        seq = seq[i + line_length:]
    if len(seq) > 0:
        print(seq)


def find_clump(sequence, interval, kmer):

    seqlen = len(sequence)
    start  = 0
    max_count = 0
    result = []

    saved_interval = interval
    while start < seqlen:
        if interval > len(sequence[start:]):
            interval = len(sequence[start:])
        count, max_strings_found = frequent_words(sequence[start:start+interval], kmer)
        #print (start, count, max_strings_found)
        if count > max_count:
            result.append([start, count])
            max_count = count
        start += interval

    '''
    interval = saved_interval
    start = int(interval / 2)
    while start < seqlen:
        if interval > len(sequence[start:]):
            interval = len(sequence[start:])
        count, max_strings_found = frequent_words(sequence[start:start + interval], kmer)
        print(start, count, max_strings_found)
        if count >= max_count:
            result.append(start, count, max_strings_found)
        start += interval
    '''

    return result



def skew_diagram(sequence):

    index = 0
    skew = 0
    skew_min_index = 0
    skew_max_index = 0
    skew_min = 0
    skew_max = 0

    for index in range(len(sequence)):
        if sequence[index] == 'C':
            skew -= 1
        elif sequence[index] == 'G':
            skew += 1
        if skew < skew_min:
            skew_min = skew
            skew_min_index = index
        elif skew > skew_max:
            skew_max = skew
            skew_max_index = index
    return skew_min, skew_min_index, skew_max, skew_max_index



def hamming_distance(seq1, seq2):

    ham = 0
    seqlen = len(seq1)
    if  seqlen != len(seq2):
        return (-1)
    else:
        i = 0
        for i in range(seqlen):
            if seq1[i] != seq2[i]:
                ham += 1
        return ham



def approxiate_match(sequence, pattern, ham):

    count = 0
    list = []
    for i in range(len(sequence)-len(pattern)):
        pattern1 = sequence[i:i+len(pattern)]
        if hamming_distance(pattern, pattern1) <= ham:
            count += 1
            list.append(pattern1)
    return count, list


def main():

    genome_vc = SeqIO.read('vibrio_cholerae.fasta.txt', 'fasta')
    ori_vc = SeqIO.read('ori_vibrio_cholerae.fasta.txt', 'fasta')


    #
    # Test the pattern_count function to determine how often a pattern occurs in a sequence
    #

    print('Start of test script\n\n')

    print_header(7,'Determine how often a pattern occurs in a sequence')

    seq = "GACCATCAAAACTGATAAACTACTTAAAAATCAGT"
    pattern = "A"
    n, pos = pattern_count(seq, pattern)
    print(f'Pattern {pattern} was found {n} times in positions {pos}.\n')

    pattern = "AA"
    n, pos = pattern_count(seq, pattern)
    print(f'Pattern {pattern} was found {n} times in positions {pos}.\n')

    pattern = "AAA"
    n, pos = pattern_count(seq, pattern)
    print(f'Pattern {pattern} was found {n} times in positions {pos}.\n')

    pattern = "AAAA"
    n, pos = pattern_count(seq, pattern)
    print(f'Pattern {pattern} was found {n} times in positions {pos}.\n')

    pattern = "AAAAA"
    n, pos = pattern_count(seq, pattern)
    print(f'Pattern {pattern} was found {n} times in positions {pos}.\n')

    pattern = "AAAAAA"
    n, pos = pattern_count(seq, pattern)
    print(f'Pattern {pattern} was found {n} times in positions {pos}.\n')

    pattern = "GATA"
    n, pos = pattern_count(seq, pattern)
    print(f'Pattern {pattern} was found {n} times in positions {pos}.\n')

    #
    # Test the frequent_words function to find the occurrence of kmers of size k
    #

    print_header(10, 'Find the occurrence of 3, 5, 7, 9 kmers of size k')

    k = 3
    count, max_strings_found = frequent_words(str(ori_vc.seq), k)
    print(f'The maximum time a {k}-kmer was found is {count}')
    print(f'The number of patterns for that maximum is {len(max_strings_found)}')
    print('The positions are:\n')
    for l in range(len(max_strings_found)):
        print('\t', max_strings_found[l])
    print('\n')

    k = 5
    count, max_strings_found = frequent_words(str(ori_vc.seq), k)
    print(f'The maximum time a {k}-kmer was found is {count}')
    print(f'The number of patterns for that maximum is {len(max_strings_found)}')
    print('The positions are:\n')
    for l in range(len(max_strings_found)):
        print('\t', max_strings_found[l])
    print('\n')

    k = 7
    count, max_strings_found = frequent_words(str(ori_vc.seq), k)
    print(f'The maximum time a {k}-kmer was found is {count}')
    print(f'The number of patterns for that maximum is {len(max_strings_found)}')
    print('The positions are:\n')
    for l in range(len(max_strings_found)):
        print('\t', max_strings_found[l])
    print('\n')

    k = 9
    count, max_strings_found = frequent_words(str(ori_vc.seq), k)
    print(f'The maximum time a {k}-kmer was found is {count}')
    print(f'The number of patterns for that maximum is {len(max_strings_found)}')
    print('The positions are:\n')
    for l in range(len(max_strings_found)):
        print('\t', max_strings_found[l])
    print('\n')

    #
    # Test the find_max_GC function to find in what block the GC percentage is highest
    #

    print_header(0, 'Find in what block the GC percentage is highest')
    seq = genome_vc.seq

    maxperc, maxindex = find_max_GC(seq)
    print(f'\nThe highest GC count in the sequence is: {maxperc:4.2f}')

    maxperc, maxindex = find_max_GC(seq, 500)
    print(f'\nThe highest GC count in the sequence cut up in {500} is: {maxperc:4.2f} at index {maxindex}\n\n')

    #
    # Test the reverse_complement function
    #

    print_header(11, 'Test the reverse_complement function')

    myseq = Seq("ATGATCAAG", IUPAC.unambiguous_dna)
    print(f'The reverse complement of {myseq} is: {myseq.reverse_complement()}\n')


    #
    # Page 13 How many times does ATGATCAAG occur in Vibrio cholerae?
    #

    print_header(13, 'How many times does ATGATCAAG occur in Vibrio cholerae?')

    pattern = 'ATGATCAAG'
    n, positions = pattern_count(genome_vc.seq,pattern)
    print(f'The patten \'{pattern}] occurs {n} times, on positions')
    print(positions)

    #
    # Page 14 Find the 9-mers in ori_thermotoga
    #


    print_header(14, 'Find the 9-mers in ori_thermotoga')
    k = 9
    ori_t = SeqIO.read('ori_thermotoga1.fasta.txt', 'fasta')
    ss = str(ori_t.seq).upper()

    count, max_strings_found = frequent_words(ss, k)
    print(f'The maximum time a {k}-kmer was found in the thermotoga ori is {count}')
    print(f'The number of patterns for that maximum is {len(max_strings_found)}')
    print('The positions are:\n')
    ll = len(max_strings_found)
    for l in range(len(max_strings_found)):
        print('\t', max_strings_found[l])
    print('\n')

    if (False):

    #
    # Page 14 Find the highest number of k-mers in blocks of length 500
    #

        interval = 1000
        print_header(14, f'Find the highest number of k-mers in blocks of length {interval}')

        sequence = str(genome_vc.seq)
        result = find_clump(sequence, interval , 8)
        max9mer = 0
        for r in result:
            if r[1] > max9mer:
                max9mer = r[1]
                index = r[0]
        count, max_strings_found = frequent_words(sequence[index:index + interval], 8)
        print(f'The maximum time a {k}-kmer was found is {count}')
        print(f'The number of patterns for that maximum is {len(max_strings_found)}')
        print(f'The start position of the block is {index}')
        print('The positions are:\n')
        for l in range(len(max_strings_found)):
            print('\t', max_strings_found[l])
        print('\n')

    print_header(0, 'Skew')

    mins, minindex, maxskew, maxindex = skew_diagram(genome_vc.seq)

    n1 = 1


    #
    # Page 27 Hamming distance
    #

    print_header(27, 'Hamming distance')

    str1 = "AACGTCCCATTC"

    str2 = "AACGTCCCATTC"
    print(f'Hamming distance of {str1} and {str2} = {hamming_distance(str1, str2)}')

    str2 = "AACGTCCCATTT"
    print(f'Hamming distance of {str1} and {str2} = {hamming_distance(str1, str2)}')

    str2 = "AACAAAACATTT"
    print(f'Hamming distance of {str1} and {str2} = {hamming_distance(str1, str2)}')

    #
    # Page 28 Approximate match
    #

    print_header(28, 'Approximate match')

    #print(approxiate_match('AACAAGCATAAACATTAAAGAG', 'AAAAA', 0))

    pattern = 'AAAAA'
    ham = 1
    n, plist = approxiate_match('AACAAGCATAAACATTAAAGAG', pattern, ham)
    print (f'There are {n} partial matches with hamming distance {ham}')
    print (f'The matches are: {plist}\n')

    pattern = 'AAAAA'
    ham = 2
    n, plist = approxiate_match('AACAAGCATAAACATTAAAGAG', pattern, ham)
    print(f'There are {n} partial matches with hamming distance {ham}')
    print(f'The matches are: {plist}\n')

    pattern = 'AAAAA'
    ham = 3
    n, plist = approxiate_match('AACAAGCATAAACATTAAAGAG', pattern, ham)
    print(f'There are {n} partial matches with hamming distance {ham}')
    print(f'The matches are: {plist}\n')

    pattern = 'AAAAA'
    ham = 4
    n, plist = approxiate_match('AACAAGCATAAACATTAAAGAG', pattern, ham)
    print(f'There are {n} partial matches with hamming distance {ham}')
    print(f'The matches are: {plist}\n')



if __name__ == "__main__":
    main()
