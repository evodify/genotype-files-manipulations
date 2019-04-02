#!/usr/bin/python2

"""
This script phases the sequences by random split of heterozygous sites.

# input:

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

# output:

CHROM   POS sample2_1   sample2_2   sample7_1   sample7_2 
chr_1   1   N   N   N   N
chr_1   2   T   C   C   C
chr_1   3   -   -   C   C
chr_1   4   T   T   T   T
chr_1   6   C   C   C   C
chr_2   1   A   A   A   A
chr_2   2   C   C   C   C
chr_2   3   N   N   N   N
chr_2   4   T   T   C   C
chr_2   5   C   C   T   T
chr_3   1   N   N   N   N
chr_3   2   C   C   C   C
chr_3   3   N   N   N   N
chr_3   4   T   T   T   T
chr_3   5   N   N   C   C
chr_4   1   N   N   C   C
chr_4   2   N   N   C   C
# command:

$ python pseudoPhasingHetero.py -i inputfile -o outputfile -s "sample2,sample7"

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


counter = 0

print('Opening the file...\n')

with open(args.input) as datafile:
  header_words = datafile.readline().split()

  ChrPos = header_words[0:2]
  ChrPosP = '\t'.join(str(e) for e in ChrPos)

  # index samples
  sampCol = calls.indexSamples(sampleNames, header_words)

  # create lists for output
  sampColnames = calls.selectSamples(sampCol, header_words)

  # create slitted column names
  idsheader = []
  for i in sampColnames:
    idsheader.append("%s_1\t%s_2 " % (i, i))
  header = '\t'.join(str(e) for e in idsheader)

  print('Creating the output file...\n')
  outname = '%s' % args.output
  outfile = open(outname, 'w')
  outfile.write('%s\t%s\n' % (ChrPosP, header))  

  print('Randomizing alleles in heterozygous genotypes ...\n')

  for line in datafile:
    words = line.split()
    chr_pos = words[0:2]

    # create lists for output
    alleles = calls.selectSamples(sampCol, words)

    # check if all genotypes are correct
    calls.if_all_gt_correct(alleles, line)

    # phase
    phasedGT = calls.pseudoPhase(alleles)
    phasedGTsplit = []
    for gt in phasedGT:
            GT = gt.split('/')
            phasedGTsplit.append(GT[0])
            phasedGTsplit.append(GT[1])

    chromPosP = '\t'.join(str(e) for e in chr_pos)
    phasedGTsplitP = '\t'.join(str(e) for e in phasedGTsplit)
    outfile.write("%s\t%s\n" % (chromPosP, phasedGTsplitP))
    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
outfile.close()
print('\nDone!')