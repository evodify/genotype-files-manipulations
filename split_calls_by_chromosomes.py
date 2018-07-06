#!/usr/bin/python2

"""
This script splits calls file into several files by chromosomes.

# input file:

CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   1   A   W   N   N   A   N   N   N   N
chr_1   2   C   Y   Y   N   C   C   N   C   N
chr_1   3   C   N   -   N   C   C   C   C   C
chr_1   4   T   T   T   N   T   T   T   T   T
chr_1   6   C   N   C   N   C   C   C   C   C
chr_2   1   A   A   A   N   A   A   A   A   A
chr_2   2   C   C   C   N   C   C   C   C   C
chr_2   3   C   N   N   N   N   N   N   N   N
chr_2   4   C   C   T   C   C   C   C   C   C
chr_2   5   T   T   C   T   Y   T   Y   T   T
chr_3   1   G   G   N   N   G   N   N   N   N
chr_3   2   C   S   C   N   C   C   N   C   N
chr_3   3   N   N   N   N   N   N   N   N   N
chr_3   4   N   T   T   N   T   T   T   T   N
chr_3   5   G   -   N   N   G   G   G   C   G
chr_4   1   G   -   N   N   G   G   G   C   G
chr_4   2   G   -   N   N   G   G   G   C   G

# test.out_chr_1.tab:
POS	POS	REF	sample1	sample2	sample3	sample4	sample5	sample6	sample7	sample8
chr_1	1	A	W	N	N	A	N	N	N	N
chr_1	2	C	Y	Y	N	C	C	N	C	N
chr_1	3	C	N	-	N	C	C	C	C	C
chr_1	4	T	T	T	N	T	T	T	T	T
chr_1	6	C	N	C	N	C	C	C	C	C

# test.out_chr_2.tab:
POS	POS	REF	sample1	sample2	sample3	sample4	sample5	sample6	sample7	sample8
chr_2	1	A	A	A	N	A	A	A	A	A
chr_2	2	C	C	C	N	C	C	C	C	C
chr_2	3	C	N	N	N	N	N	N	N	N
chr_2	4	C	C	T	C	C	C	C	C	C
chr_2	5	T	T	C	T	Y	T	Y	T	T

# test.out_chr_3.tab:
POS	POS	REF	sample1	sample2	sample3	sample4	sample5	sample6	sample7	sample8
chr_3	1	G	G	N	N	G	N	N	N	N
chr_3	2	C	S	C	N	C	C	N	C	N
chr_3	3	N	N	N	N	N	N	N	N	N
chr_3	4	N	T	T	N	T	T	T	T	N
chr_3	5	G	-	N	N	G	G	G	C	G

# test.out_chr_4.tab:
POS	POS	REF	sample1	sample2	sample3	sample4	sample5	sample6	sample7	sample8
chr_4	1	G	-	N	N	G	G	G	C	G
chr_4	2	G	-	N	N	G	G	G	C	G



# command:

$ python split_calls_by_chromosomes.py -i test.tab -o test.out

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process (optional)', type=str, required=False)
args = parser.parse_args()

# check if samples names are given and if all sample names are present in a header
sampleNames = calls.checkSampleNames(args.samples, args.input)

############################# program #############################

CHRprev = 'NA'

with open(args.input) as datafile:
    header_words = datafile.readline().split()

    # index a sample
    sampCol = calls.indexSamples(sampleNames, header_words)
    header_samples = calls.selectSamples(sampCol, header_words)
    header_samplesP = '\t'.join(str(e) for e in header_samples)

    for line in datafile:
        words = line.split()
        CHR = words[0]
        Pos = words[1]

        samples = calls.selectSamples(sampCol, words)
        samplesP = '\t'.join(str(e) for e in samples)

        # find chromosome border
        if CHRprev == CHR:
            output.write('%s\t%s\t%s\n' % (CHR, Pos, samplesP))
        else:
            CHRprev = CHR
            print "Writing", CHR, "to the file ..."
            output = open(args.output + "_" + str(CHR) + ".tab", 'w')
            output.write('%s\t%s\t%s\n' % (header_words[1], header_words[1], header_samplesP))
            output.write('%s\t%s\t%s\n' % (CHR, Pos, samplesP))

datafile.close()
output.close()
