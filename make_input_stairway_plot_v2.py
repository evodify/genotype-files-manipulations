#! /usr/bin/env python
'''
This script extracts the information required for a [stairway plot v2](https://sites.google.com/site/jpopgen/stairway-plot) input.

input.file:

#CHROM xPOS sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8 sample9 sample10 sample11 sample12 sample13 sample14 sample15 sample16 sample17 sample18 sample19 sample20 sample21 sample22
scaffold_1 585 C C C A C N N C C C C C C C C C C C C C C C
scaffold_1 586 G G G G G G G G G G G G G G G G G G G G G G
scaffold_1 587 G G C C G G G G G G G G G G G G G G G G G G
scaffold_1 589 C C C C C C C C C C C C C C C C C C C C C C
scaffold_1 591 G G G G G G G G G G G G G G G G G G G G G G
scaffold_1 593 T T T T T T A A A A A A A A A A A A A A A A
scaffold_1 594 T T T T G G G G G G G G G G G G G G G G G G
scaffold_1 595 G G G G G G G G G G G G G G G G G G G G G G
scaffold_1 596 C C C C C C C C C C C C C C C C C C C C C C
scaffold_1 597 T T T T T T T T T T T T T T T T T T T T T T

ancestor.file:

#CHROM  xPOS    Ancestor
scaffold_1  585 A
scaffold_1  586 G
scaffold_1  587 G
scaffold_1  589 C
scaffold_1  591 G
scaffold_1  593 T
scaffold_1  594 T
scaffold_1  595 G
scaffold_1  596 T
scaffold_1  597 A

Output:

nseq:	22
L:	9
SFS:	1	1	1

These lines need to be placed in an input blueprint file for the stairway plot v2.

# command

$ python make_input_stairway_plot_v2.py -i input.file -o output.file -a ancestor.file

Optionally samples can be specified with -s "sample1,sample2...sampleN"

To use sites with missing data, specify number of allowed Ns. This will produced reduced SFS by using subsampling from hypergeometric distribution. (Hernandez, Ryan D., et al. "Context-dependent mutation rates may cause spurious signatures of a fixation bias favoring higher GC-content in humans." Molecular biology and evolution 24.10 (2007): 2196-2202.)

contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu
'''

############################# modules #############################

import calls # my custom module
import collections
import numpy as np # to use missing data

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-a', '--ancestor', help = 'name of the file with ancestral sequence to polarize alleles', type=str, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process (optional)', type=str, required=False)
parser.add_argument('-m', '--missing', help = 'allowed number of missing data', type=int, required=False)
args = parser.parse_args()

# check if samples names are given and if all sample names are present in a header
sampleNames = calls.checkSampleNames(args.samples, args.input)

# check the missing data threshold
if args.missing:
  AlowedN = int(args.missing)
else:
  AlowedN = 0
  print 'The option "-m" is not specified. All sites with missing data will be skipped.'

############################# program #############################

ref = open(args.ancestor, 'r')
ref_header = ref.readline()
words2 = ref.readline().split()
ref_chr_pos = words2[0:2]
ref_ch = int(ref_chr_pos[0].split('_')[1])
ref_pos = int(ref_chr_pos[1])
ancest = words2[2]

output = open(args.output, 'w')

seqL = 0
counter = 0
freq = []

with open(args.input) as datafile:
  header_line = datafile.readline()
  header_words = header_line.split()

  # index samples
  sampCol = calls.indexSamples(sampleNames, header_words)
  nseq = len(sampleNames)
  maxGenot = nseq - AlowedN
  output.write("nseq:\t%s\n" % maxGenot)

  for line in datafile:
    words = line.split()
    ch = int(words[0].split('_')[1])
    pos = int(words[1])

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

    # select samples
    genotypesN = calls.selectSamples(sampCol, words)

    # skip sites with missing data
    numN =  calls.countPerPosition(genotypesN, 'N')
    if numN <= AlowedN: 
      genotypes = [i for i in genotypesN if i != 'N']
      nGenotypes = len(genotypes)
    else:
      continue # skip lines with too many Ns 

    # count alleles
    numAl = collections.Counter(genotypes)
    numAlM = numAl.most_common()

    while ch > ref_ch or (ch == ref_ch and pos > ref_pos):
      words2 = ref.readline().split()
      if words2 == []:
        ancest = 'N'
        break
      else:
        ref_chr_pos = words2[0:2]
        ref_ch = int(ref_chr_pos[0].split('_')[1])
        ref_pos = int(ref_chr_pos[1])
        ancest = words2[2]

    if ancest == 'N':  # skip unpolarized sites
      continue
    elif len(numAl) == 1:  # if fixed
      al1 = numAlM[0][0]
      if al1 == ancest:
        freq.append(0)
      else:
        n1 = numAlM[0][1]
        x1 = np.random.hypergeometric(n1, nGenotypes-n1, maxGenot, 1)[0] # subsample
        freq.append(x1)
    elif len(numAl) == 2:  # if biallelic
      al1 = numAlM[0][0]
      al2 = numAlM[1][0]
      if al2 == ancest:
        if nGenotypes == maxGenot:
          x1 = numAlM[0][1]
        else:
          n1 = numAlM[0][1]
          x1 = np.random.hypergeometric(n1, nGenotypes-n1, maxGenot, 1)[0] # subsample
        freq.append(x1)
      elif al1 == ancest:
        if nGenotypes == maxGenot:
          x2 = numAlM[1][1]
        else:
          n2 = numAlM[1][1]
          x2 = np.random.hypergeometric(n2, nGenotypes-n2, maxGenot, 1)[0] # subsample
        freq.append(x2)
    seqL+=1

  output.write("L:\t%s\n" % seqL)
  SFSnum = collections.Counter(freq)
  SFSnumS = collections.OrderedDict(sorted(SFSnum.items()))
  output.write("SFS:")
  print "\nSFS:"
  for key, value in SFSnumS.items():
    print key, value
    if key == 0 or key == maxGenot:
      continue
    else:
      output.write("\t%s" % (value))

datafile.close()
output.close()
ref.close()
print('Done!')
