#! /usr/bin/python
'''
This script assesses a number of SNPs with a sliding window approach.

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


# output.tab:

CHR	midPOS	numberSNPs
chr1	5	6
chr1	15	1
chr2	5	5
chr2	15	1

# command

$ python numberSNPsPerWindows.py \
-i input.tab \
-w 10 \
-o output

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

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
    '-w',
    '--window',
    help = 'window size',
    type=int,
    required=True)
args = parser.parse_args()

############################# program #############################

def processWindow(CHR, POS, numberSNPs, output):
    output.write('%s\t%s\t%s\n' % (CHR, POS, numberSNPs))

output = open(args.output, 'w')
output.write('CHR\tmidPOS\tnumberSNPs\n')

lineNumber = 0
CHRprev = 'NA'
POSprev = 'NA'
numberSNPs = 'NA'
windows_step = args.window
windowsEnd = windows_step
halfWindow = float(windows_step)/2

with open(args.input) as datafile:
    datafile.readline()
    for line in datafile:
        words = line.split()
        CHRcurrent = words[0]
        POScurrent = int(words[1])

        if CHRprev != CHRcurrent:
            if numberSNPs != 'NA':
                POS = int(windowsEnd - halfWindow)
                processWindow(CHRprev, POS, numberSNPs, output)
            numberSNPs = 1
            windowsEnd = windows_step
        elif POScurrent > windowsEnd:
            POS = int(windowsEnd - halfWindow)
            processWindow(CHRprev, POS, numberSNPs, output)
            numberSNPs = 1
            windowsEnd += windows_step
        else:
            numberSNPs += 1
        CHRprev = CHRcurrent
        POSprev = POScurrent

        lineNumber = calls.lineCounter(lineNumber)
    POS = int(windowsEnd - halfWindow)
    processWindow(CHRprev, POS, numberSNPs, output)

print 'Done!'
