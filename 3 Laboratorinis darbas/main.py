import math
from collections import Counter
from bioinfokit.analys import fastq
from Bio import SeqIO
from Bio.Blast import NCBIWWW
	
RANGES = {
    'Sanger': (33, 73),
    'Illumina-1.8': (33, 74),
    'Solexa': (59, 104),
    'Illumina-1.3': (64, 104),
    'Illumina-1.5': (67, 104)
}


def checkEncoding(encoding, encoding_name):
    value = True
    i = len(unique_ascii_symbols)
    while i > 0:
        if unique_ascii_symbols[i - 1] not in encoding:
            print("Is not encoded in ", encoding_name, " because contains character: ", unique_ascii_symbols[i - 1])
            value = False
            break
        i -= 1
    return value


def listToString(s):
    str1 = ""
    for ele in s:
        str1 += ele
    return str1


def GCratio(sequence):
    letterC = sequence.count('C')
    letterG = sequence.count('G')
    ratio = (letterC + letterG) / float(len(sequence))
    return round(ratio, 2)


def round(number, decimals=0):
    multiplier = 10 ** decimals
    return math.ceil(number * multiplier) / multiplier


if __name__ == "__main__":
    qual_list = []
    GCRatio_list = []
	fastq_iter = fastq.fastq_reader(file='reads_for_analysis.fastq')
    
    for record in fastq_iter:
        header_1, sequence, header_2, qual = record
		
        # sequence_len = len(sequence)
        # a_base = sequence.count('A')
        # print(sequence, qual, a_base, sequence_len)
        
		qual_list.append(qual)
        # print(header_1, sequence, header_2, qual)
        
		GCRatio_list.append(GCratio(sequence))
        
		if 0.24 <= GCratio(sequence) <= 0.39:
            print("First peak sequence: ", sequence)
			
        if 0.47 <= GCratio(sequence) <= 0.58:
            print("Second peak sequence: ", sequence)
			
        if 0.64 <= GCratio(sequence) <= 0.76:
            print("Third peak sequence: ", sequence)

    unique_ascii_symbols = list(set(listToString(qual_list)))
    print("Unique ASCII symbols: ", unique_ascii_symbols)

    sanger = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI"""
    solexa = ";<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"
    illumina_1_3 = "@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"
    illumina_1_5 = "CDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghi"
    illumina_1_8 = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"""

    if checkEncoding(sanger, "sanger"):
        print("The encoding is done in sanger")
    if checkEncoding(solexa, "solexa"):
        print("The encoding is done in solexa")
    if checkEncoding(illumina_1_3, "illumina 1.3+"):
        print("The encoding is done in illumina 1.3+")
    if checkEncoding(illumina_1_5, "illumina 1.5+"):
        print("The encoding is done in illumina 1.5+")
    if checkEncoding(illumina_1_8, "illumina 1.8+"):
        print("The encoding is done in illumina 1.8+")

    counter = Counter(GCRatio_list)
	
    print("All CGRatio's")
    print(counter)

	print("End of the program...")
