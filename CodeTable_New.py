import pprint
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq, MutableSeq


class CodeTable:

    GeneticCode = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
                   "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
                   "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
                   "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
                   "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
                   "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                   "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
                   "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
                   "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
                   "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                   "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
                   "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
                   "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
                   "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                   "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
                   "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

    FromThreeToOneTable = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                           'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                           'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                           'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    FromOneToThreeTable = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
                           'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN',
                           'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP',
                           'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

    TranslationCountTable = {}
    ReverseGeneticCode = {}
    TranslationCountTable = {}
    ReverseGeneticCode = {}

    def __init__(self):

        TranslationCountTable = {}
        ReverseGeneticCode = {}
        TranslationCountTable = {}
        ReverseGeneticCode = {}

        for codon in self.GeneticCode:
            am_acid = self.GeneticCode[codon]
            if am_acid not in self.TranslationCountTable:
                self.ReverseGeneticCode[am_acid] = [codon]
                self.TranslationCountTable[am_acid] = 1
            else:
                self.TranslationCountTable[am_acid] += 1
                self.ReverseGeneticCode[am_acid].append(codon)


    def three_letter_code(self, am1):
        if am1 in self.FromOneToThreeTable:
            return self.FromOneToThreeTable[am1]
        else:
            return ''


    def one_letter_code(self, am3):
        if am3 in self.FromThreeToOneTable:
            return self.FromThreeToOneTable[am3]
        else:
            return ''


    def get_protein(self, seq):
        aminostring = ""
        alen = len(seq)
        alen = alen // 3

        try:
            for i in range(alen):
                str1 = seq[i * 3:i * 3 + 3]
                aa = self.GeneticCode[str1]
                if aa != 'STOP':
                    aminostring += aa
                else:
                    return aminostring
        except:
            pass
            return ''

        return aminostring

    def get_dna(self, aminoseq):

        codon_list = []
        i = 0
        rna_sequences = []

        # First get the codons that code for each amino acids
        for aa in aminoseq:
            codons = self.ReverseGeneticCode[aa.upper()]
            codon_list.append(self.ReverseGeneticCode[aa.upper()])
            i += 1

        # Then make the combinations of these codons
        nr_aminoacids = len(codon_list)
        combine_list = []
        base_list = codon_list[0]
        if nr_aminoacids == 1:
            base_list = codon_list[0]
        else:
            for i in range(1, nr_aminoacids):
                for bl in base_list:
                    for al in codon_list[i]:
                        new = bl + al
                        combine_list.append(new)
                base_list = combine_list[:]
                combine_list = []

        rna_list = base_list

        # Now change the found RNA sequences into DNA sequences
        dna_sequences = []
        for rna in rna_list:
            str1 = rna.replace('U', 'T')
            dna_sequences.append(str1)
        return dna_sequences


def how_many_combinations(seq):
    am_acids = seq.split('-')
    n = 1
    for aa in am_acids:
        n = n * TranslationCountTable[aa.upper()]
    return n


def add_next_level(base_list, list_to_add):
    combine_list = []

    for bl in base_list:
        for al in list_to_add:
            new = bl + al
            combine_list.append(new)

    return combine_list


def find_coding_sequences(dna_sequence, amino_sequence):
    codon_list = []
    i = 0
    rna_sequences = []

    # First get the codons that code for each amino acids
    for aa in amino_sequence:
        codons = ReverseGeneticCode[aa.upper()]
        codon_list.append(ReverseGeneticCode[aa.upper()])
        i += 1

    # Then make the combinations of these codons
    nr_aminoacids = len(codon_list)
    rna_list = codon_list[0]
    for i in range(1, nr_aminoacids):
        rna_list = add_next_level(rna_list, codon_list[i])

    # Now change the found RNA sequences into DNA sequences
    dna_sequences = []
    for rna in rna_list:
        str1 = rna.replace('U', 'T')
        dna_sequences.append(str1)

    # Then see if you can find it in the given DNA sequence
    myseq = Seq(dna_sequence, IUPAC.unambiguous_dna)
    for dna in dna_sequences:
        if myseq.find(dna) != -1:
            print(f'Found {dna}')

    # Then see if you can find it in the reverse complement DNA sequence
    myseq = myseq.reverse_complement()
    for dna in dna_sequences:
        if myseq.find(dna) != -1:
            print(f'Found {dna}')

def print_long_string(sequence, width):
    seqlen = len(sequence)

    index = 0
    while seqlen > 0:
        prlen = min(width, seqlen)
        print(sequence[index:index + prlen])
        index += prlen
        seqlen -= prlen


if __name__ == "__main__":

    ct = CodeTable()
    s1 = ct.one_letter_code("VAL")
    s2 = ct.one_letter_code("SER")
    s3 = ct.three_letter_code("Y")
    s4 = ct.three_letter_code("W")
    s5 = ct.get_protein('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA')
    s5 = ct.get_protein('GAAACU')
    print(s5)
    s5 = ct.get_protein('AGUUUC')
    print(s5)

    s7 = ct.get_dna('MA')
    myseq = []
    for dna in s7:
        myseq.append(Seq(dna))
        myseq.append(Seq(dna).reverse_complement())


    print(len(s7))

    pp = pprint.PrettyPrinter(width=10, compact=True)
    pp.pprint (s7)
    n=1

