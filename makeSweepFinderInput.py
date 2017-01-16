#!/usr/bin/python2

"""
This script makes an input file for [SweepFinder](http://people.binf.ku.dk/rasmus/webpage/sf.html). To make such file information on the outgroup/ancestral sequence is required.

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
scaffold_1  117 T
scaffold_1  122 C
scaffold_1  137 A
scaffold_1  139 T
scaffold_1  148 T
scaffold_1  161 C
scaffold_1  170 C
scaffold_1  174 A

# fai:

chr_1      19624517        12      80      81
chr_2      14106692        19869848        80      81
chr_3      15060676        34152886        80      81

# output:

position    x   n   folded
117 3   3   0
137 1   3   0
148 2   3   0
170 1   3   0

# command:

$ python makeSweepFinderInput.py -i inputfile -o outputfile -f file.fai -a ancestral.tab -N 1

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""
############################# modules #############################

import calls # my custom module
import collections
import re
import random

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-a', '--ancestor', help = 'name of the outgroup/ancestral sequence file', type=str, required=True)
parser.add_argument('-f', '--fai', help = 'name of the fasta.fai file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-N', '--missing', help = 'number of allowed Ns', type=int, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process (optional)', type=str, required=False)
args = parser.parse_args()

# check if samples names are given and if all sample names are present in a header
sampleNames = calls.checkSampleNames(args.samples, args.input)

############################# program #############################

counter = 0

Ns = args.missing
fai = open(args.fai, 'r')
fai_words = fai.readline().split()
fai_ch = fai_words[0]
fai_start1 = 0
fai_start2 = int(fai_words[1])

output = open(args.output, 'w')
output.write("position\tx\tn\tfolded\n")
outputChr = open(fai_ch+'_'+args.output, 'w')
outputChr.write("position\tx\tn\tfolded\n")

ancestFile = open(args.ancestor, 'r')
ancest_header = ancestFile.readline()
words2 = ancestFile.readline().split()
ancest_chr_pos = words2[0:2]
ancest_ch = int(ancest_chr_pos[0].split('_')[1])
ancest_pos = int(ancest_chr_pos[1])
ancest = words2[2]

print('Opening the file...')
with open(args.input) as datafile:
  header_line = datafile.readline()
  header_words = header_line.split()

  # index samples
  sampCol = calls.indexSamples(sampleNames, header_words)

  for line in datafile:
    words = line.split()
    chr_pos = words[0:2]
    ch = int(words[0].split('_')[1])
    pos = int(words[1])

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

    # select alleles
    alleles = calls.selectSamples(sampCol, words)

    # count missing data
    numAlN = collections.Counter(alleles)
    valueN = numAlN['N']

    if valueN <= Ns:  # filer by missing data threshold
      AllallesNoN = [i for i in alleles if i != 'N']
    else:
      continue

    # count alleles
    numAlNoN = collections.Counter(AllallesNoN)
    numAl = numAlNoN.most_common()
 
    # find overlap with ancestral sequence
    while (ch > ancest_ch) or (ch == ancest_ch and pos > ancest_pos):
      words2 = ancestFile.readline().split()
      if words2 == []:
        ancest = 'N'
        break
      else:
        ancest_chr_pos = words2[0:2]
        ancest_ch = int(ancest_chr_pos[0].split('_')[1])
        ancest_pos = int(ancest_chr_pos[1])
        ancest = words2[2]

    # find overlap with fai file to define chromosome borders
    if ch <= 8:  # major chromosomes that will be split
      while chr_pos[0] != fai_ch:
        fai_words = fai.readline().split()
        if fai_words == []:
          break
        else:
          fai_ch = fai_words[0]
          fai_start1 = fai_start1+fai_start2
          fai_start2 = int(fai_words[1])
          outputChr = open(fai_ch+'_'+args.output, 'w')
          outputChr.write("position\tx\tn\tfolded\n")
      posP = pos+fai_start1
    #  concatenate small scaffolds
    elif ch == 9:
      while chr_pos[0] != fai_ch:
        fai_words = fai.readline().split()
        if fai_words == []:
          break
        else:
          fai_ch = fai_words[0]
          fai_start1 = fai_start1+fai_start2
          fai_start2 = int(fai_words[1])
          outputChr = open("scaffolds_small"+'_'+args.output, 'w')
          outputChr.write("position\tx\tn\tfolded\n")
      posP = pos+fai_start1
    else:
      while chr_pos[0] != fai_ch:
        fai_words = fai.readline().split()
        if fai_words == []:
          break
        else:
          fai_ch = fai_words[0]
          fai_start1 = fai_start1+fai_start2
          fai_start2 = int(fai_words[1])
      posP = pos+fai_start1

    # polarize alleles
    n = len(AllallesNoN)
    al1 = numAl[0][0]
    x1 = numAl[0][1]
    f = 0

    if len(numAl) > 2: # skip biallelic
      continue
    elif chr_pos == ancest_chr_pos:  # folded
      if len(numAl) == 1: # fixed
        if ancest == 'N':
          f = 1
          x = random.choice([x1, 0])
        elif al1 == ancest:
          x = 0
        else:
          x = n
      elif len(numAl) == 2: # biallelic
        al2 = numAl[1][0]
        x2 = numAl[1][1]
        if al1 == ancest:
          x = x2
        elif al2 == ancest:
          x = x1
        else:
          f = 1
          x = random.choice([x1, x2])
    else:    # unfolded
      if len(numAl) == 1: # fixed
        f = 1
        x = random.choice([x1, 0])
      elif len(numAl) == 2:
        x2 = numAl[1][1]
        f = 1
        x = random.choice([x1, x2])

    if x == 0 and f == 0:  # skip sites with fixed ancestral alleles
      continue
    else:
      output.write("%s\t%s\t%s\t%s\n" % (posP, x, n, f))
      outputChr.write("%s\t%s\t%s\t%s\n" % (pos, x, n, f))

datafile.close()
output.close()
outputChr.close()
fai.close()
ancestFile.close()
print('Done!')
