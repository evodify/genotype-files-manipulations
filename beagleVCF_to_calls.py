#! /usr/bin/python
'''
This script converts Beagle phased VCF to calls file.

# beagle.vcf

#CHROM  POS   ID   REF  ALT QUAL FILTER  INFO FORMAT sample1  sample2
chr6    215    .    T    G    .    PASS    .    GT    1|0    1|0
chr6    245    .    A    C    .    PASS    .    GT    0|0    0|0
chr6    261    .    C    A    .    PASS    .    GT    0|0    0|0
chr6    280    .    C    A    .    PASS    .    GT    0|0    1|0
chr6    321    .    G    C    .    PASS    .    GT    0|0    0|0
chr6    328    .    C    A    .    PASS    .    GT    0|0    0|0
chr6    345    .    G    A    .    PASS    .    GT    0|0    0|0
chr6    365    .    A    G    .    PASS    .    GT    0|0    0|0
chr6    367    .    G    C    .    PASS    .    GT    0|0    0|0
chr7    1215    .    T    G    .    PASS    .    GT    1|0    1|0
chr7    1245    .    A    C    .    PASS    .    GT    0|0    0|0
chr7    1261    .    C    A    .    PASS    .    GT    0|0    0|0
chr7    1280    .    C    A    .    PASS    .    GT    0|0    1|0
chr7    1321    .    G    C    .    PASS    .    GT    0|0    0|0
chr7    1328    .    C    A    .    PASS    .    GT    0|0    0|0
chr7    1345    .    G    A    .    PASS    .    GT    0|0    0|0
chr7    1365    .    A    G    .    PASS    .    GT    0|0    0|0
chr7    1367    .    G    C    .    PASS    .    GT    0|0    0|0

# output.tab

CHROM  POS	sample1	sample2
chr6	215	G|T	G|T
chr6	245	A|A	A|A
chr6	261	C|C	C|C
chr6	280	C|C	A|C
chr6	321	G|G	G|G
chr6	328	C|C	C|C
chr6	345	G|G	G|G
chr6	365	A|A	A|A
chr6	367	G|G	G|G
chr7	1215	G|T	G|T
chr7	1245	A|A	A|A
chr7	1261	C|C	C|C
chr7	1280	C|C	A|C
chr7	1321	G|G	G|G
chr7	1328	C|C	C|C
chr7	1345	G|G	G|G
chr7	1365	A|A	A|A
chr7	1367	G|G	G|G


# command:

$ python beagleVCF_to_calls.py \
-i beagle.vcf \
-o output.tab

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls  # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument(
    '-i',
    '--input',
    help='name of the input file',
    type=str,
    required=True)
parser.add_argument(
    '-o',
    '--output',
    help='name of the output file',
    type=str,
    required=True)
args = parser.parse_args()

############################# program #############################

output = open(args.output, 'w')

counter = 0

with open(args.input) as datafile:

    header = datafile.readline().split()
    samplesNames = header[9:]
    samplesNamesP = '\t'.join(str(s) for s in samplesNames)  
    output.write("CHROM\tPOS\t%s\n" % samplesNamesP)

    for line in datafile:
        words = line.split()
        chr = words[0]
        pos = int(words[1])
        ref = words[3]
        alt = words[4]
        gt = words[9:]

        genotypes = []
        for g in gt:
            gtSplit = g.split("|")
            gt1 = calls.convert01toATGC(ref, alt, gtSplit[0])
            gt2 = calls.convert01toATGC(ref, alt, gtSplit[1])
            gt12 = '|'.join(gt1+gt2) 
            genotypes.append(gt12)

        genotypesP = '\t'.join(str(g) for g in genotypes)  
        output.write("%s\t%s\t%s\n" % (chr, pos, genotypesP))

        counter += 1
        if counter % 100000 == 0:
            print str(counter), "lines processed"

datafile.close()
output.close()
print('Done!')
