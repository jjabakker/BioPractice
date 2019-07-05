


#
#
#
#
#
#
#
#

def str_overlap(str1, str2):

    sol1 = str_overlap_work(str1, str2)
    sol2 = str_overlap_work(str2, str1)

    if sol1[0] == sol2[0] and sol1[0] != 0:
        # There are two similar overlaps
        return sol1[0], sol1[1], sol2[1]
    if sol1[0] > sol2[0]:
        # There is an overlap with str1 before str2
        return sol1[0], '', sol1[1]
    elif sol1[0] < sol2[0]:
        # There is an overlap with str2 before str1
        return sol2[0], '', sol2[1]
    else:
        # There is no overlap
        return sol1


def str_overlap_work(str1, str2):

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



def main():
    n = str_overlap('abcde', 'cabc')
    n = str_overlap('abcdef', 'bcdefghi')
    n = str_overlap('abcdef', 'defghi')
    n = str_overlap('abcdef', 'efghi')
    n = str_overlap('abcdef', 'fghi')
    n = str_overlap('abcdef', 'wbcd')

    n = str_overlap('abcde', 'cccccabc')
    n = str_overlap('abcdeccc', 'ccccabc')



    n = str_overlap('abcdeccc', 'bcd')
    n = str_overlap('abcdeccc', 'bcdde')
    n = str_overlap('abcdeccc', 'bcdec')
    n = str_overlap('bcdec', 'abcdeccc')

    n = str_overlap('AAACGTAACT', 'ACTTTAAGG')

    n = 1
if __name__ == "__main__":
    main()
