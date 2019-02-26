#! /usr/bin/python
'''
This script converts genotype calls file to .geno and .snps files 
suitable for XPCLR (https://reich.hms.harvard.edu/software)
# input.tab
CHROM	POS	ANC	DER	sample1	sample2	sample3	sample4	sample5	sample6	sample7
chr1	1	A	T	0/1	./.	./.	0/0	./.	./.	./.
chr1	2	C	T	1/0	0/1	./.	0/0	0/0	1/1	0/0
chr1	3	C	-	./.	1/1	./.	0/0	0/0	0/0	0/0
chr1	4	G	T	1/1	1/1	./.	1/1	1/1	1/1	1/1
chr1	6	C	T	./.	0/0	./.	0/0	0/0	0/0	0/0
chr2	1	A	T	0/0	0/0	./.	0/0	0/0	0/0	0/0
chr2	2	T	C	1/1	1/1	./.	1/1	1/1	1/1	1/1
chr2	4	C	T	0/0	1/1	0/0	0/0	0/0	0/0	0/0
chr2	5	T	C	0/0	1/1	0/0	1/1	0/0	1/0	0/0
chr3	1	G	C	0/0	./.	./.	0/0	./.	./.	./.
chr3	2	C	G	1/0	0/0	./.	0/0	0/0	./.	0/0
chr3	4	C	T	1/1	1/1	1/1	1/1	1/1	1/1	1/1
You can get this input file from your call file with polarize_callsTab.py
You can also use "|" separator fro phased data.
# family1output.geno:
0 1 9 9
1 0 0 1
9 9 1 1
0 0 1 1
0 0 1 1
1 0 0 0
# family2output.geno:
9 9 9 9
0 0 1 1
0 0 0 0
0 0 0 0
0 0 1 0
0 0 9 9
# output.snps:
chr1_1 chr1 1e-06 1 A T
chr1_2 chr1 2e-06 2 C T
chr1_3 chr1 3e-06 3 C -
chr2_4 chr2 4e-06 4 C T
chr2_5 chr2 5e-06 5 T C
chr3_2 chr3 2e-06 2 C G
NOTE: By default I assume 1 bp == 1e-06 cM.
You can adjust the cM column with make_custom_map.py
# command
$ python calls_to_XPCLR_input.py \
-i input.tab \
-o output \
-f "family1[sample1,sample2];family2[sample5,sample6]"
# contact:
Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu
'''

############################# modules #############################

import calls # my custom module
import re # to split input
import random # for randomization

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument(
  '-i',
  '--input',
  help = 'name of the input file',
  type = str,
  required = True)
parser.add_argument(
  '-o',
  '--output',
  help = 'name of the output file',
  type = str,
  required = True)
parser.add_argument(
  '-f',
  '--family',
  help = 'Specify the family list in the format \
          "family1[sample1,sample2];family2[sample5,sample6]"',
  type=str,
  required=True)
args = parser.parse_args()

print "Checking samples name ..."
familyNames = args.family
Fsamples = []
familySamples = {}
for i in familyNames.strip("\"").split(";"):
  famName = i.split("[")[0]
  famSample = re.split("\[|\]|", i)[1]
  Fsamples.append(famSample.split(","))
  familySamples[famName] = calls.checkSampleNames(famSample, args.input)
samples = calls.flattenList(Fsamples)
calls.checkSampleNames('ANC', args.input)
calls.checkSampleNames('DER', args.input)

############################# program #############################

outputSNPs = open(args.output + '.snps', 'w')

with open(args.input) as datafile:
    header_words = datafile.readline().split()
    sampleIndex = calls.indexSamples(samples, header_words)
    ANCindex = calls.indexSamples(['ANC'], header_words)
    DERindex = calls.indexSamples(['DER'], header_words)

    FamilyIndex = {}
    phasedLines = {}
    for family in familySamples:
        FamilyIndex[family] = calls.indexSamples(familySamples[family],
                                                 header_words)
        phasedLines[family] = []

    for line in datafile:
        words = line.split()
        CHR = words[0]
        POS = words[1]
        chr_pos = CHR + '_' + POS
        ANC = calls.selectSamples(ANCindex, words)[0]
        DER = calls.selectSamples(DERindex, words)[0]
        allSeq = calls.selectSamples(sampleIndex, words)

        if calls.is_biallelic(allSeq):
            for family in FamilyIndex:
                famSeq = calls.selectSamples(FamilyIndex[family],  words)
                GTsplit = []
                splitfamSeq = []
                for GT in famSeq:
                    GT = GT.replace('/', ' ')
                    GT = GT.replace('|', ' ')
                    GT = GT.replace('.', '9')
                    GTsplit.append(GT)
                splitSeq = ' '.join(str(e) for e in GTsplit)
                splitfamSeq.append(splitSeq)
                phasedLines[family].append(splitfamSeq)
            cM = float(POS)*0.000001
            outputSNPs.write('%s %s %s %s %s %s\n' %
                            (chr_pos, CHR, cM, POS, ANC, DER))
        else:
            print('WARNING: Site %s is not polymorphic or biallelic. '
                    'Skipping it ...' % (chr_pos))

    for family in familySamples:
        outputGeno = open(family + '_' + args.output + '.geno', 'w')
        for line in phasedLines[family]:
            lineP = ' '.join(str(e) for e in line)
            outputGeno.write("%s\n" % lineP)
        outputGeno.close()  

print 'Done!'
