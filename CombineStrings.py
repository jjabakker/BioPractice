


#
# The function tries to get a maximum string length
#
#     str_overlap('abcabc', 'abcabc') gives abcabcabc
#
#     abcaaa and aaagfh: abcaaagfh or abcaaaagfh or abcaaaaagfh ?
#
#

def str_overlap(str1, str2):

    sol1 = str_overlap_work(str1, str2)
    sol2 = str_overlap_work(str2, str1)

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


def str_overlap_work_old(str1, str2):

    n = len(str1)

    i = 0

    for i in range(n):

        s1 = str1[i:]
        s2 = str2

        n = min(len(s1),len(s2))

        if (s1[:n] == str2[:n]):
            index_in_str1 = i
            s3 = str1 + str2[n:]
            return n, s3

    return  0, ''

def str_overlap_work(str1, str2):

    n = len(str1)
    nmax  = 0
    smax = ''

    i = 0

    for i in range(n-1,-1,-1):

        s1 = str1[i:]
        s2 = str2

        n = min(len(s1),len(s2))
        s11 = s1[:n]
        s22 = str2[:n]
        if (s1[:n] == str2[:n]):
            index_in_str1 = i
            s3 = str1 + str2[n:]
            if n > nmax :
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


def main():
    s1 = 'aaabbb'
    s2 = 'aaaccc'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)

    s2 = 'cccaaa'
    n = str_overlap(s1, s2)  # OK. No overlap
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

    s1 = 'abcabc'
    s2 = 'abcabc'
    n = str_overlap(s1, s2)  # OK. No overlap
    print_res(n, s1, s2)

    

    n = str_overlap('AAACGTAACT', 'ACTTTAAGG')

    n = 1
if __name__ == "__main__":
    main()
