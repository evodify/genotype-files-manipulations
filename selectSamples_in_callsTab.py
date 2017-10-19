#!/usr/bin/env python2

"""
This script subsamples a genotype calls file by sample names. It also can be used to rearrange samples in a calls file.

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

# output:

CHROM   POS   sample3   sample4   sample6   sample8   REF
chr_1   1   N   A   N   N   A
chr_1   2   N   C   N   N   C
chr_1   3   N   C   C   C   C
chr_1   4   N   T   T   T   T
chr_1   6   N   C   C   C   C
chr_2   1   N   A   A   A   A
chr_2   2   N   C   C   C   C
chr_2   3   N   N   N   N   C
chr_2   4   C   C   C   C   C
chr_2   5   T   Y   Y   T   T
chr_3   1   N   G   N   N   G
chr_3   2   N   C   N   N   C
chr_3   3   N   N   N   N   N
chr_3   4   N   T   T   N   N
chr_3   5   N   G   G   G   G
chr_4   1   N   G   G   G   G
chr_4   2   N   G   G   G   G


# command:

$ python2 selectSamples.py -i input.tab -o output -s "sample3,sample4,sample6,sample8,REF"

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

print('Opening the file...')

with open(args.input) as datafile:
  header_words = datafile.readline().split()

  # index samples
  sampCol = calls.indexSamples(sampleNames, header_words)
  
  # make output header
  print('Creating the output file...')
  fileoutput = open(args.output, 'w')
  sampHeader = calls.selectSamples([0,1]+sampCol, header_words)
  sampHeaderP = '\t'.join(str(el) for el in sampHeader)
  fileoutput.write(sampHeaderP+'\n')

  for line in datafile:
    words = line.split()
    chr_pos = words[0:2]

    # select samples
    genotypes = calls.selectSamples(sampCol, words)

    # make output
    chr_posP = '\t'.join(str(el) for el in chr_pos)
    genotypesP = '\t'.join(str(el) for el in genotypes)
    fileoutput.write('%s\t%s\n' % (chr_posP, genotypesP))

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
fileoutput.close()

print('Done!')
