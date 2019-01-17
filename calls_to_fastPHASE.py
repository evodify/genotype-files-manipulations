#!/usr/bin/python2
"""
This script converts genotype calls file to PHASE/fastPHASE input.
Input should contain information only for one chromosome.

PHASE/fastPHASE: http://stephenslab.uchicago.edu/software.html

!!! 
I experienced problems with data sets larger than 300K SNPs.
I think it was a problem of fastPHASE, because smaller dataset converted 
with this script worked fine.
!!!

# input file:

CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   1   A   W   N   A   A   A   T   T   T
chr_1   2   C   Y   Y   N   C   C   T   C   T
chr_1   3   C   -   -   C   C   C   C   C   C
chr_1   4   T   T   T   N   T   T   T   T   T
chr_1   6   C   N   C   N   C   C   C   C   C

# output:

9
5
P 1 2 3 4 6
SSSSS
#REF
ACCTC
ACCTC
#sample1
TT?T?
AC?T?
#sample2
?T?TC
?C?TC
#sample3
A?C??
A?C??
#sample4
ACCTC
ACCTC
#sample5
ACCTC
ACCTC
#sample6
TTCTC
TTCTC
#sample7
TCCTC
TCCTC
#sample8
TTCTC
TTCTC

# command:

$ python2 calls_to_fastPHASE.py \
    -i input.tab \
    -o output \
    -s "REF,sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8"

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

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
parser.add_argument(
    '-s',
    '--samples',
    help='column names of the samples to process (optional)',
    type=str,
    required=False)
args = parser.parse_args()

sampleNames = calls.checkSampleNames(args.samples, args.input)

############################# program #############################

callsDF = calls.callsParser(args.input, sampleNames)

numChromosomes = set(callsDF.chrmosomes)
if len(numChromosomes) != 1:
  raise IOError('Input should contain information only for one chromosome. '
                'You input has %s chromosomes. '
                'Please, split the input by chromosomes.' %
                len(numChromosomes))

output = open(args.output, 'w')

NumberPos = len(callsDF.positions)
NumberSamp = len(sampleNames)

# write  no.individuals and no.SNPsites
output.write('%s\n%s\nP' % (NumberSamp, NumberPos))

prevPos = 0
for pos in callsDF.positions:
    if int(prevPos) > int(pos):
        raise IOError('Positions should be increasing. '
                      'Positions %s goes after %s' % (pos, prevPos))
    prevPos = pos

posChunk = calls.chunks(callsDF.positions, 500000)
sssChunk = calls.chunks(callsDF.positions, 500000)

# write postions
for posCh in posChunk:
    posChP = ' '.join(str(i) for i in posCh)
    output.write(" %s" % posChP)
output.write("\n")

# indicate that all postions are SNPs
for sssCh in sssChunk:
    SSS = ['S']*len(posCh)
    sssP = ''.join(str(i) for i in SSS)
    output.write("%s" % sssP)
output.write("\n")


for s in callsDF.names:
    # write ID
    output.write("#%s\n" % s)
    # write genotypes
    haplo1 = []
    haplo2 = []
    for g in  callsDF[s]:
        gTwo = calls.pseudoPhase(g)[0].split('/')
        for i in range(len(gTwo)):
            if gTwo[i] in 'N-.':
                gTwo[i] = '?'

        haplo1.append(gTwo[0])
        haplo2.append(gTwo[1])
    
    haplo1c = calls.chunks(haplo1, 500000)
    haplo2c = calls.chunks(haplo2, 500000)

    for chunk1 in haplo1c:
        haplo1P = ''.join(str(i) for i in chunk1)
        output.write("%s\n" % haplo1P)

    for chunk2 in haplo2c:
        haplo2P = ''.join(str(i) for i in chunk2)
        output.write("%s\n" % haplo2P)

output.close()
