from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.Data import CodonTable
from Bio import SeqIO
from Bio import Align
import Bio


#
#
#

def test_seq():
    print(f'\n\n{Bio.__version__}')

    myseq = Seq("GATCGAAATGGGCCTAAAAATATAGGATCGAAAATCGC", IUPAC.unambiguous_dna)

    print(myseq.alphabet)
    print(myseq)
    print(myseq.__len__())

    # You can calculate the GC ratio yourself or you can you use a special function to do so
    print(f'The number of G\'s is: {myseq.count("G")}')
    print(f'The number of C\'s is: {myseq.count("C")}')
    print(100 * (myseq.count("G") + myseq.count("C")) / myseq.__len__())
    print(f'The ratio GC in the string is: {GC(myseq)}')

    # Sequences can be inverted and complemented
    print(f'The string is:\t\t\t\t {myseq}')
    print(f'The reverse is:\t\t\t\t {myseq[::-1]}')
    print(f'The complement is:\t\t\t {myseq.complement()}')
    print(f'The reverse complement is:\t {myseq.reverse_complement()}')

    # A simple function to determine if a string is palindromic
    print(f'GAAG is palindrome: {ispalindromic("GAAG")}')
    print(f'GAAG is palindrome: {ispalindromic("GTAG")}')

    # Find all the sequences of AAA
    # Start with the first AAA and that start looping
    positions = []
    pos = myseq.find('AAA')
    while pos != -1:
        positions.append(pos)
        pos = myseq.find('AAA', pos + 1)
    print(positions)
    count = len(positions)

    # A demo to show that Seq is not mutable.
    # If you need to change a Seq, make it mutable first
    try:
        myseq[0] = 'C'

    except:
        print("Oops! As Seq is immutable ")

    im_seq = myseq.tomutable()
    im_seq[0] = 'C'

    seq = im_seq
    print(seq)


def test_translation():
    """

    Both RNA and DNA can be translated....
    """

    coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", IUPAC.unambiguous_dna)
    template_dna = coding_dna.reverse_complement()
    messenger_RNA = coding_dna.transcribe()

    print('\n\n')
    print(f'The coding DNA is    \t\t\t\t\t\t: {coding_dna}')
    print(f'The template DNA is  \t\t\t\t\t\t: {template_dna} (Remember it is reversed and displayed as 5-3')
    print(
        f'The messenger RNA is \t\t\t\t\t\t: {messenger_RNA} (The mRNA is the coding string with the T replaced by U)')
    print(f'Now if you translate the mRNA\t\t\t\t: {messenger_RNA.translate(table="Standard")}')
    print(f'Same result if you translate the coding DNA\t: {coding_dna.translate()}')
    print(f'You can stop at first stop codon\t\t\t: {messenger_RNA.translate(to_stop=True)}')

    # For the translation you can specify which tables you want to use
    print('\n\n')
    print(CodonTable.unambiguous_dna_by_name['Standard'])
    print('\n\n')
    print(CodonTable.unambiguous_dna_by_name['Vertebrate Mitochondrial'])


def test_simple_file_io():
    """
    Here read a file with single record in one go
    """

    record = SeqIO.read('NC_005816.fna.txt', 'fasta')
    print(record)
    print(f'The record name is:        {record.name}')
    print(f'The record description is: {record.description}')
    print(f'The record sequence is:    {record.seq}')
    pass


def test_file_io_fasta():
    for seq_record in SeqIO.parse('ls_orchid.fasta.txt', 'fasta'):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))

    for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
        print(seq_record.seq)

    records = list(SeqIO.parse('ls_orchid.fasta.txt', 'fasta'))

    print('\n\nThe first record\n')
    print(records[0])

    print('\n\nThe second record\n')
    print(records[1])


def test_file_io_gbk():
    for seq_record in SeqIO.parse('ls_orchid.gbk.txt', 'genbank'):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))

    for seq_record in SeqIO.parse('ls_orchid.gbk.txt', 'genbank'):
        print(seq_record + '\n\n')


def ispalindromic(sequence):
    return sequence == sequence[::-1]


def test_align():
    """
    Determining how a string matches is a key challenge...

    """

    aligner = Align.PairwiseAligner()
    searchstr = 'ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG'
    searchfor = 'GAATGG'
    alignments = aligner.align(searchstr, searchfor)
    for alignment in alignments:
        print(alignment)

    score = aligner.score
    print(score)
    print(aligner)


def main():
    test_seq()
    test_translation()
    test_simple_file_io()
    test_file_io_fasta()
    test_file_io_gbk()
    test_align()


if __name__ == "__main__":
    main()
