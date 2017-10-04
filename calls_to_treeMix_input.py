#!/usr/bin/env python2
'''
This script outputs alleles count file that is required as input for TreeMix https://bitbucket.org/nygcresearch/treemix/wiki/Home.

inputfile:

CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   1   A   A   A   A   A   A   A   N   T
chr_1   2   N   C   C   T   C   C   T   C   C
chr_1   3   T   A   A   T   T   C   T   A   T
chr_1   4   N   N   A   C   C   C   A   C   C


outputfile:

pop1    pop2
4,0 2,1
3,1 3,1
1,2 1,3

command:

python calls_to_treeMix_input.py -i test.input -o test.out -p "pop1[sample1,sample2,sample3,sample4];pop2[sample5,sample6,sample7,sample8]"

contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu
'''

############################# modules #############################

import calls # my custom module
import re
from  collections import Counter

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-p', '--pop', help = 'Specify the populations in the format "pop1[sample1,sample2,];pop2[sample5,sample6]"', type=str, required=True)
#parser.add_argument('-p', '--parental', help = 'Specify a parental groups', type=str, required=True)

args = parser.parse_args()

# check and append population names and samples
popNames = args.pop
pops = []
for popi in popNames.strip("\"").split(";"):
  popName = popi.split("[")[0]
  popSample = re.split("\[|\]", popi)[1]
  pops.append(popName)
  vars()[popName + "samples"] = calls.checkSampleNames(popSample,args.input)

############################## program #############################

counter = 0

print('Opening the file...')

with open(args.input) as datafile:
  header_words = datafile.readline().split()

  # make output header
  print('Creating the output file...')
  fileoutput = open(args.output, 'w')
  popsP = '\t'.join(str(w) for w in pops)
  fileoutput.write("%s\n" % popsP)

  for popName in pops:
    # index samples
    vars()[popName + "Index"]  = calls.indexSamples(vars()[popName + "samples"], header_words)

  for line in datafile:
    words = line.split()
    GT = words[2:]
    GTpair = [i for i in list(set(GT)) if i != 'N'] # get the set of alleles, skip missing alleles
    popNum = 0

    #if ("N" not in GT) and (len(GTpair) == 2) : # skip missing data or non-biallelic
    if (len(GTpair) == 2): # skip non-biallelic
      for popName in pops:
        popNum += 1 # to make correct output. See below
        # select genotypes
        sGT = calls.selectSamples(vars()[popName + "Index"], words)

        counts = Counter(sGT) # count alleles
        bicounts = [counts[GTpair[0]],counts[GTpair[1]]] # extract counts
        bicountsP = ','.join(str(w) for w in bicounts)

        # make output
        if popNum == len(pops):
          fileoutput.write("%s\n" % bicountsP)
        else:
          fileoutput.write("%s\t" % bicountsP)

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
fileoutput.close()

print('Done!')
