#! /usr/bin/python
'''
This script converts genotype calls file to .geno and .snps files 
suitable for XPCLR (https://reich.hms.harvard.edu/software)

# input.tab

CHROM	POS	ANC	DER	sample1	sample2	sample3	sample4	sample5	sample6	sample7
chr1	1	A	T	0/1	./.	./.	0/0	./.	./.	./.
chr1	2	C	T	1/0	0/1	./.	0/0	0/0	1/1	0/0
chr1	3	C	-	./.	1/1	./.	0/0	0/0	0/0	0/0
chr1	4	G	T	0/1	1/1	1/1	1/1	1/1	1/1	1/1
chr1	6	C	T	./.	0/0	./.	0/0	0/0	0/0	0/0
chr1	10	A	T	0/0	0/0	1/1	0/0	0/0	0/0	0/0
chr1	12	A	T	0/0	0/0	0/1	0/0	0/0	0/0	0/0
chr2	1	T	C	1/1	1/1	1/1	0/1	1/1	1/1	1/1
chr2	4	C	T	0/0	1/1	0/0	0/0	0/0	0/0	0/0
chr2	5	T	C	0/0	1/1	0/0	1/1	0/0	1/0	0/0
chr2	8	G	C	0/0	1/1	1/1	0/0	1/1	1/1	1/1
chr2	9	C	G	1/0	0/0	./.	0/0	0/0	./.	0/0
chr2	12	C	T	0/0	1/1	1/1	1/1	1/1	1/1	1/1


You can get this input file from your call file with polarize_callsTab.py
You can also use "|" separator for phased data.


# family1output.geno:

0 1 9 9 9 9
1 0 0 1 9 9
9 9 1 1 9 9
0 0 0 0 1 1
0 0 0 0 0 1
1 1 1 1 1 1
0 0 1 1 0 0
0 0 1 1 0 0
0 0 1 1 1 1
1 0 0 0 9 9
0 0 1 1 1 1


# family2output.geno:

0 0 9 9 9 9
0 0 0 0 1 1
0 0 0 0 0 0
0 0 0 0 0 0
0 0 0 0 0 0
0 1 1 1 1 1
0 0 0 0 0 0
1 1 0 0 1 0
0 0 1 1 1 1
0 0 0 0 9 9
1 1 1 1 1 1


# output.snps:

chr1_1 chr1 0.0 1 A T
chr1_2 chr1 2e-06 2 C T
chr1_3 chr1 3e-06 3 C -
chr1_10 chr1 1e-05 10 A T
chr1_12 chr1 1.2e-05 12 A T
chr2_1 chr2 0.0 1 T C
chr2_4 chr2 4e-06 4 C T
chr2_5 chr2 5e-06 5 T C
chr2_8 chr2 8e-06 8 G C
chr2_9 chr2 9e-06 9 C G
chr2_12 chr2 1.2e-05 12 C T


NOTE: By default I assume 1 bp == 1e-06 cM.
You can adjust the cM column with make_custom_map.py

# command

$ python calls_to_XPCLR_input.py \
-i input.tab \
-f "family1[sample1,sample2,sample3];family2[sample4,sample5,sample6]" \
-w 5 -n 3 \
-o output

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module
import re # to split input
from random import sample # for randomization

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
parser.add_argument(
    '-w',
    '--window',
    help = 'window size',
    type=int,
    required=False)
parser.add_argument(
    '-n',
    '--numberSNPsInWindow',
    help = 'allowed number of SNPs in a window',
    type=int,
    required=False)
args = parser.parse_args()

print "Checking sample names ..."
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

############################# functions #################################

def createWindow(info, seq, wInfo, wSeq):
    if calls.is_biallelic(seq):
        windowInfoR = [info]
        windowSeqR = [seq]
    else:
        windowInfoR = []
        windowSeqR = []
        # print('WARNING: Site %s:%s is neither polymorphic nor biallelic. '
        #     'Skipping it ...' % (info[0], info[1]))
    return windowInfoR, windowSeqR

def processWindow(windowInfoI, windowSeqI, windowsSize,
                  familySamplesI, FamilyIndexI, outputSNPs):
    
    phasedLines = {}
    for family in familySamplesI:
        phasedLines[family] = []
        
    numberSNPs = len(windowInfoI)

    if windowsSize == 'NA':
        subsetIndex = range(numberSNPs)
    elif numberSNPs > windowsSize:
        subsetIndex = sorted(sample(xrange(numberSNPs), windowsSize))
    else:
        subsetIndex = range(numberSNPs)
    
    for i in subsetIndex:
        info = windowInfoI[i]
        seq = windowSeqI[i]
        CHR = info[0]
        POS = info[1]
        chr_pos = CHR + '_' + POS
        ANC = info[2]
        DER = info[3]

        for family in familySamplesI:
            famSeq = calls.selectSamples(FamilyIndexI[family],  seq)
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
        if int(POS) == 1:
            cM = 0.0
        else:
            cM = float(POS)*0.000001
        outputSNPs.write('%s %s %s %s %s %s\n' %
                        (chr_pos, CHR, cM, POS, ANC, DER))

    for family in familySamples:
        outputGeno = open(family + '_' + args.output + '.geno', 'a')
        for pline in phasedLines[family]:
            plineP = ' '.join(str(e) for e in pline)
            outputGeno.write("%s\n" % plineP)
        outputGeno.close() 

############################# program #############################

outputSNPs = open(args.output + '.snps', 'w')

if args.window and args.numberSNPsInWindow:
    windows_step = args.window
    windows_size = args.numberSNPsInWindow
else:
    windows_step = 100000
    windows_size = 'NA'

with open(args.input) as datafile:
    header_words = datafile.readline().split()
    sampleIndex = calls.indexSamples(samples, header_words)

    FamilyIndex = {}
    for family in familySamples:
        FamilyIndex[family] = calls.indexSamples(familySamples[family],
                                                 header_words[4:])
        outputGeno = open(family + '_' + args.output + '.geno', 'w')
        outputGeno.close() 
    
    windowsEnd = windows_step
    windowInfo = []
    windowSeq = []
    CHRprev = 'NA'
    POS = ''

    lineNumber = 0

    for line in datafile:
        words = line.split()
        CHRcurrent = words[0]
        POScurrent = int(words[1])
        infoLine = words[0:4]
        seqLine = calls.selectSamples(sampleIndex, words)

        if CHRprev != CHRcurrent and windowInfo:
            processWindow(windowInfo, windowSeq, windows_size,
                          familySamples, FamilyIndex, outputSNPs)
            windowInfo, windowSeq  = createWindow(infoLine, seqLine,
                                                  windowInfo, windowSeq)
            windowsEnd = windows_step
        elif POScurrent > windowsEnd:
            processWindow(windowInfo, windowSeq, windows_size,
                          familySamples, FamilyIndex, outputSNPs)
            windowInfo, windowSeq  = createWindow(infoLine, seqLine,
                                                  windowInfo, windowSeq)
            windowsEnd += windows_step
        else:
            if calls.is_biallelic(seqLine):
                windowInfo.append(infoLine)
                windowSeq.append(seqLine)
            # else:
            #     print('WARNING: Site %s:%s is neither polymorphic nor'
            #           ' biallelic. Skipping it ...' % (CHRcurrent, POScurrent))
        CHRprev = CHRcurrent

        lineNumber = calls.lineCounter(lineNumber)

    processWindow(windowInfo, windowSeq, windows_size,
                          familySamples, FamilyIndex, outputSNPs)

print 'Done!'
