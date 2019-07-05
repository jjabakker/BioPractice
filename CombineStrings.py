import random
import copy



#
# The function tries to get a maximum string length
#
# If there are multiple overlaps, the longest overlap will be chosen, not necessarily leading to the
# longest combined string
#
# The return consists of four parametrs
#
#       Int: The number cof overlap solutions (1, 2 or 3)
#       Int: The length of the overlap
#       Str: Overlap solution 1
#       Str: Overlap solution 2


def str_overlap(str1, str2, min_overlap=3):
    sol1 = str_overlap_work(str1, str2, min_overlap)
    sol2 = str_overlap_work(str2, str1, min_overlap)

    if sol1[0] == sol2[0] and sol1[0] != 0:
        # There are two similar overlaps
        return 2, sol1[0], sol1[1], sol2[1]
    elif sol1[0] > sol2[0]:
        # There is an overlap with str1 before str2
        return 1, sol1[0], sol1[1], ''
    elif sol1[0] < sol2[0]:
        # There is an overlap with str2 before str1
        return 1, sol2[0], sol2[1], ''
    else:
        # There is no overlap
        return 0, 0, '', ''


def str_overlap_work(str1, str2, min_overlap):
    n = len(str1)
    nmax = 0
    smax = ''

    for i in range(n - 1, -1, -1):

        s1 = str1[i:]
        s2 = str2

        n = min(len(s1), len(s2))
        if n >= min_overlap:
            if s1[:n] == str2[:n]:
                s3 = str1 + str2[n:]
                if n > nmax:
                    smax = s3
                    nmax = n

    return nmax, smax


def print_res(n, s1, s2):
    if n[0] == 0:
        print(f'There is no overlap in {s1} and {s2}\n')
    elif n[[0] == 1]:
        print(f'Strings {s1} and {s2} combine to {n[2]} with an overlap of {n[1]}\n')
    else:
        print(f'Strings {s1} and {s2} combine to {n[2]} and {n[3]} with an overlap of {n[1]}\n')



def construct_poem(poem, n_seg, seg_len, min_overlap):


    poem_len = len(poem)

    segments = {}
    segments_copy = {}
    del_list = {}
    not_solved = True
    n_seg = 100

    for i in range(n_seg):
        start =  random.randint(1, poem_len-15)
        seglen = random.randint(int(0.5*seg_len), seg_len)
        seg_end = min(start+seglen, poem_len)
        segments[i] = poem[start:seg_end]
        print (segments[i])

    for i in range(10):
        if i != 0:
            segments_copy = {}
            for s in segments:
                if segments[s] != '':
                    segments_copy[s] = segments[s]
            segments = copy.deepcopy(segments_copy)
        if len(segments_copy) == 1:
            return segments_copy

        for s1 in segments:
            for s2 in segments:
                if s1 != s2:
                    if segments[s1] == segments[s2]:
                        segments[s1] = ''
                    elif segments[s1] == '' or segments[s2] == '':
                        pass
                    else:
                        n = str_overlap(segments[s1], segments[s2],min_overlap)
                        if n[0] > 0 and n[1] > 2:
                            segments[s1] = ''
                            segments[s2] = n[2]
                            n_seg += 1

    return segments


def print_long_string(sequence, width):
    seqlen = len(sequence)

    index = 0
    while seqlen > 0:
        prlen = min(width, seqlen)
        print(sequence[index:index + prlen])
        index += prlen
        seqlen -= prlen


def main():
    print('\n')
    s1 = 'aaabbb'
    s2 = 'aaaccc'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)

    s2 = 'cccaaa'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)

    s2 = 'cccaaa'
    n = str_overlap(s1, s2, 5)  # OK. No overlap
    print_res(n, s1, s2)

    s2 = 'bbbaaa'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)

    s2 = 'aaccc'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)

    s2 = 'aabbb'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)

    s2 = 'cccaa'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)

    s2 = 'bbbaa'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)

    s2 = 'bbbaa'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)

    s1 = 'abcde'
    s2 = 'abcde'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)

    s1 = 'abcdeccc'
    s2 = 'ccccabc'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)

    s1 = 'abcdeccc'
    s2 = 'bcd'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)

    s2 = 'bcdec'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)

    s1 = 'bcdec'
    s2 = 'abcdeccc'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)


    s1 = 'AAACGTAACT'
    s2 = 'ACTTTAAGG'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)

    poem1 = 'We do not know what purpose — if any — these other 9-mers serve in the E. coli genome, but we do know '
    poem2 = 'that there are many different types of hidden messages in genomes; these hidden messages have a tendency '
    poem3 = 'to cluster within a genome, and most of them have nothing to do with replication.'
    poem = poem1 + poem2 + poem3

    print_long_string(poem, 80)
    print('\n')
    n = construct_poem(poem, 30, 40, 7)


    print ('\n')
    for l in n:
        print_long_string (n[l],80)




if __name__ == "__main__":
    main()
