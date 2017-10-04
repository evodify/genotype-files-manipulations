#! /usr/bin/env python

"""
This script makes input for MSMC (https://github.com/stschiff/msmc).

# input:

#CHROM  POS sample1 sample2 sample3 sample4
scaffold_1  113 G   G   G   N
scaffold_1  117 C   C   N   C
scaffold_1  122 C   C   N   C
scaffold_1  137 A   A   T   N
scaffold_1  139 T   T   T   N
scaffold_1  148 A   A   T   N
scaffold_1  161 C   C   C   T
scaffold_1  170 C   T   C   N
scaffold_1  174 A   A   A   N 

# output:

1 137 4 AAT?
1 148 2 AAT?
1 161 1 CCCT
1 170 1 CTC?

# command:

$ python make_input_MSMC.py -i datafile -o outputfile -s "sample1,sample2,sample3,sample4"


contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

############################# modules #############################

#import collections
import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-m', '--missing', help = 'missing data threshold to remove sites', type=int, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process (optional)', type=str, required=False)

args = parser.parse_args()

# check if samples names are given and if all sample names are present in a header
sampleNames = calls.checkSampleNames(args.samples, args.input)

############################# program #############################

counter = 0
siteNumber = 1
chrPrev = str(1)

print('Opening the file...')

with open(args.input) as datafile:
  header_line = datafile.readline()
  header_words = header_line.split()

  # index samples
  sampCol = calls.indexSamples(sampleNames, header_words)

  print('Creating the output file...')
  fileoutput = open(chrPrev+'_'+args.output, 'w')

  for line in datafile:
    words = line.split()
    chr = str(words[0].split('_')[1])
    pos = words[1]

    # split chromosomes into separate files
    if chr != chrPrev:
      fileoutput.close()
      fileoutput = open(chr+'_'+args.output, 'w')
      chrPrev = chr
      siteNumber = 1

    # select samples
    genotypes = calls.selectSamples(sampCol, words)

    # count Ns
    valueN = calls.countPerPosition(genotypes, 'N')

    if valueN <= args.missing:
      genotypesMerged = ''.join(str(e) for e in genotypes)
      genotypesMergedP = genotypesMerged.replace('N', '?')
    else:
      continue

    # count the number of called sites
    if not calls.is_polymorphic(genotypes):
      siteNumber += 1
    else:
      fileoutput.write("%s\t%s\t%s\t%s\n" % (chr, pos, siteNumber, genotypesMergedP))
      siteNumber = 1

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
fileoutput.close()

print('Done!')
