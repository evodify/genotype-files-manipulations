#!/usr/bin/env python2
'''
This script outputs only unique allele of one population relative to another.

I use this approach to remove ancestral variation (named as parental below) from a population of interest.

inputfile:

CHROM   POS sample1 sample2 sample3 sample4 
chr_1   1   A   N   A   A
chr_1   2   C   C   C   T
chr_1   3   C   T   N   N

outputfile:

CHROM   POS sample3 sample4
chr_1   2   N   T

command:

python find_popSpesificAlleles.py -i test.input -o test.out -s "sample3,sample4" -p "sample1,sample2"

contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu
'''

############################# modules #############################

import calls # my custom module
import re

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s', '--sample', help = 'Specify the sample group', type=str, required=True)
parser.add_argument('-p', '--parental', help = 'Specify a parental groups', type=str, required=True)

args = parser.parse_args()

# check if samples names are given and if all sample names are present in a header
sNames = calls.checkSampleNames(args.sample, args.input)
pNames = calls.checkSampleNames(args.parental, args.input)

############################# program #############################

counter = 0

print('Opening the file...')

with open(args.input) as datafile:
  header_words = datafile.readline().split()

  # index samples
  sIndex = calls.indexSamples(sNames, header_words)
  pIndex = calls.indexSamples(pNames, header_words)

  # create lists for output
  sNames = calls.selectSamples(sIndex, header_words)

  # make output header
  print('Creating the output file...')
  fileoutput = open(args.output, 'w')
  sampGroups = '\t'.join(str(w) for w in header_words[0:2]+sNames)
  fileoutput.write("%s\n" % sampGroups)

  for line in datafile:
    words = line.split()
    chr_pos = words[0:2]

    # select samples
    sGT = calls.selectSamples(sIndex, words)
    pGT = calls.selectSamples(pIndex, words)

    # compare genotypes
    pGTset = set(pGT)
    if "N" in pGTset: # remove N from parental group
      pGTset.remove("N")
    sGPprint = []
    for sgt in sGT:
      if sgt not in pGTset and sgt != "N" and pGTset:
        sGPprint.append(sgt)
      else:
        sGPprint.append("N")

    # make output
    if not calls.all_missing(sGPprint):
      chr_posP = '\t'.join(str(el) for el in chr_pos)
      sGPprintP = '\t'.join(str(el) for el in sGPprint)
      fileoutput.write('%s\t%s\n' % (chr_posP, sGPprintP))

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
fileoutput.close()

print('Done!')
