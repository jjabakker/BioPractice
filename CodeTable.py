import pprint

GeneticCode = {"UUU": "F", "UUC": "F", "UUA": "L",    "UUG": "L",
               "UCU": "S", "UCC": "S", "UCA": "S",    "UCG": "S",
               "UAU": "Y", "UAC": "Y", "UAA": "STOP", "UAG": "STOP",
               "UGU": "C", "UGC": "C", "UGA": "STOP", "UGG": "W",
               "CUU": "L", "CUC": "L", "CUA": "L",    "CUG": "L",
               "CCU": "P", "CCC": "P", "CCA": "P",    "CCG": "P",
               "CAU": "H", "CAC": "H", "CAA": "Q",    "CAG": "Q",
               "CGU": "R", "CGC": "R", "CGA": "R",    "CGG": "R",
               "AUU": "I", "AUC": "I", "AUA": "I",    "AUG": "M",
               "ACU": "T", "ACC": "T", "ACA": "T",    "ACG": "T",
               "AAU": "N", "AAC": "N", "AAA": "K",    "AAG": "K",
               "AGU": "S", "AGC": "S", "AGA": "R",    "AGG": "R",
               "GUU": "V", "GUC": "V", "GUA": "V",    "GUG": "V",
               "GCU": "A", "GCC": "A", "GCA": "A",    "GCG": "A",
               "GAU": "D", "GAC": "D", "GAA": "E",    "GAG": "E",
               "GGU": "G", "GGC": "G", "GGA": "G",    "GGG": "G"}


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


def get_3_letter_code(am1):
    if am1 in FromOneToThreeTable:
        return FromOneToThreeTable[am1]
    else:
        return ''


def get_1_letter_code(am3):
    if am3 in FromThreeToOneTable:
        return FromThreeToOneTable[am3]
    else:
        return ''


def make_reverse_code_table():

    for codon in GeneticCode:
        am_acid = GeneticCode[codon]
        if am_acid not in TranslationCountTable:
            ReverseGeneticCode[am_acid] = [codon]
            TranslationCountTable[am_acid] = 1
        else:
            TranslationCountTable[am_acid] += 1
            ReverseGeneticCode[am_acid].append(codon)

    pp = pprint.PrettyPrinter(width=60, compact=True)
    pp.pprint(TranslationCountTable)
    pp.pprint(ReverseGeneticCode)


def translate_protein(sequence):


    aminostring = ""
    alen = len(sequence)
    alen = alen // 3

    try:
        for i in range(alen):
            str = sequence[i*3:i*3+3]
            aa = GeneticCode[str]
            aminostring = aminostring + aa
    except:
        return ''

    return aminostring


def how_cany_combinations(sequence):

    am_acids = sequence.split('-')
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

    # First get the codons for each amino acids

    for aa in amino_sequence:
        codons = ReverseGeneticCode[aa.upper()]
        codon_list.append(ReverseGeneticCode[aa.upper()])
        i += 1

    # Then make the combinations
    nr_aminoacids = len(codon_list)


    dna_list = codon_list[0]

    for i in range(1, nr_aminoacids):
        dna_list = add_next_level(dna_list, codon_list[i])

    for s1 in codon_list[0]:
        for s2 in codon_list[1]:
            str = s1 + s2
            rna_sequences.append(str)

    dna_sequences = []
    for rna in rna_sequences:
        str = rna.replace('U', 'T')
        dna_sequences.append(str)

    for dna in dna_sequences:
        if dna_sequence.find(dna):
            print(f'Found {dna}')


def print_long_string(sequence, width):

    seqlen = len(sequence)

    index = 0
    while seqlen > 0 :
        prlen = min(width, seqlen)
        print(sequence[index:index+prlen])
        index += prlen
        seqlen -= prlen



if __name__ == "__main__":
    make_reverse_code_table()

    f = open("sequence_to translate.txt", "r")
    sequence = f.read()

    print(translate_protein('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'))
    # output should be MAMAPRTEINSTRING
    print_long_string((translate_protein(sequence)),80)

    print(f'There are {how_cany_combinations("V-K-L-F-P-Y-F-N-Q-W")}, combinations for string "V-K-L-F-P-Y-F-N-Q-W"')

    m = 1

    s = get_3_letter_code('L')
    s = get_1_letter_code(s)
    s = get_3_letter_code('J')

    find_coding_sequences('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA','LALA')

    print(translate_protein('AUGGCC'))
    print(translate_protein('GGCCAU'))
    print(translate_protein('AUGGCC'))
    print(translate_protein('WWW'))

    print(translate_protein('AUGWWW'))

    print_long_string('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBBBBBBBBBBBCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCDDD',30)