#! /usr/bin/env python

"""
This script polarizes the genotype data by keeping only derived alleles relative to an outgroup/ancestral sequence.

Note! Chromosome number in the first column must be separated by _.
For example, chr_1 - correct, chr1 - incorrect.

# input:

#CHROM  POS sample1 sample2 sample3 sample4
scaffold_1  113 G   G   G   N
scaffold_1  117 C   C   N   C
scaffold_1  122 C   C   N   C
scaffold_1  137 A   A   T   N
scaffold_1  139 T   T   T   N
scaffold_1  148 A   A   T   N
scaffold_1  161 C   C   C   N
scaffold_1  170 C   T   C   N
scaffold_1  174 A   A   A   N

# ancestral:

#CHROM  POS Co_ancestor
scaffold_1  113 G
scaffold_1  117 N
scaffold_1  122 C
scaffold_1  137 A
scaffold_1  139 T
scaffold_1  148 T
scaffold_1  161 C
scaffold_1  170 C
scaffold_1  174 A

# output:

#CHROM  POS sample1 sample2 sample3 sample4
scaffold_1  137 N   N   T   N
scaffold_1  148 A   A   N   N
scaffold_1  170 N   T   N   N


# command:

$ python polarizeGT.py -i input.tab -a ancestral.tab -o output

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""
############################# modules #############################

import calls # my custom module

############################# options #############################


parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the reference input file', type=str, required=True)
parser.add_argument('-a', '--ancestral', help = 'name of the outgroup/ancestral sequence file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process (optional)', type=str, required=False)
args = parser.parse_args()

# check if samples names are given and if all sample names are present in a header
sampleNames = calls.checkSampleNames(args.samples, args.input)

############################# program #############################

counter = 0

# read the header
ances = open(args.ancestral, 'r')
ances_words = ances.readline()

output = open(args.output, 'w')

print('Opening the file...')
with open(args.input) as datafile:
  header_line = datafile.readline()
  header_words = header_line.split()

  # index samples
  sampCol = calls.indexSamples(sampleNames, header_words)

  # make a header of the output file
  print('Creating the output file...')
  samples_head = calls.selectSamples(sampCol, header_words)
  samples_headP = '\t'.join(str(e) for e in header_words[0:2]+samples_head)
  output.write('%s\n' % samples_headP)

  # read the second line of the ancestral file
  ances_words = ances.readline().split()
  ances_ch = int(ances_words[0].split('_')[1])
  ances_pos = int(ances_words[1])
  ances_gt = ances_words[2]

  for line in datafile:
    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

    words = line.split()
    ch = int(words[0].split('_')[1])
    pos = int(words[1])

    # find overlap
    while (ch > ances_ch) or (ch == ances_ch and pos > ances_pos):
      ances_words = ances.readline().split()
      if ances_words == []:
        break
      else:
        ances_ch = int(ances_words[0].split('_')[1])
        ances_pos = int(ances_words[1])
        ances_gt = ances_words[2]

    # introduce missing data if there is no overlap
    if pos != ances_pos:
      continue  # skip all missing data lines
    else:
      # select samples
      samples_gt = calls.selectSamples(sampCol, words)
      for i in range(len(samples_gt)):
        if samples_gt[i] == ances_gt or ances_gt == "N":  # replace ancestral alleles with Ns, keep derived, or skip if ancestral state is N.
          samples_gt[i] = 'N'
      if calls.all_missing(samples_gt): # skip all missing data lines
        continue
      else:
        wordsP = '\t'.join(str(e) for e in words[0:2]+samples_gt)

    # write output
    output.write('%s\n' % wordsP)



datafile.close()
ances.close()
output.close()

print('Done!')
