#!/usr/bin/python2

"""
This script cuts genotype calls file with the given window size and outputs FASTA files for every window.

# input file:

CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   1   A   W   N   A   A   N   N   N   N
chr_1   2   C   Y   Y   A   C   C   A   C   C
chr_1   3   C   N   -   A   C   C   C   C   C
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

# output1.fasta:

>sample5
NCCTC
>sample4
ACCTC
>sample7
NCCTC
>sample6
NACTC
>sample1
WYNTN
>sample3
AAANN
>sample2
NY-TC
>REF
ACCTC
>sample8
NCCTC

# output2.fasta:

>sample5
ACNCT
>sample4
ACNCY
>sample7
ACNCT
>sample6
ACNCY
>sample1
ACNCT
>sample3
NNNCT
>sample2
ACNTC
>REF
ACCCT
>sample8
ACNCT

# command:

$ python2 slidingWindowSNPs.py -i test.tab -o output -w 5 -N 3


# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

############################# modules #############################

import calls # my custom module
from collections import defaultdict

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process (optional)', type=str, required=False)
parser.add_argument('-w', '--window', help = 'window size', type=int, required=True)
parser.add_argument('-N', '--missingness', help = 'allowed number of Ns', type=int, required=True)

args = parser.parse_args()

# check if samples names are given and if all sample names are present in a header
sampleNames = calls.checkSampleNames(args.samples, args.input)
Wsize = args.window
allowedN = args.missingness

############################# program #############################

windows_count = 0
windows_size = 0
CHRprev = 'NA'

outputNames = open(args.output+'.Names', 'w')

with open(args.input) as datafile:
  header_words = datafile.readline().split()

  # index a sample
  sampCol = calls.indexSamples(sampleNames, header_words)

  for line in datafile:
    words = line.split()

    # find chromosome border
    if CHRprev != words[0]:
      genotype = []
      gend = defaultdict(list)
      windows_size = 0
      genN = [0]*len(sampleNames)
    CHRprev = words[0]
    windows_size += 1

    # append genotypes to samples dict
    G = calls.selectSamples(sampCol, words)
    for s,g,i in zip(sampleNames,G,range(len(genN))):
      gend[s].append(g)
      if g == "N":
        genN[i] += 1

    # if the window size is reached:
    if windows_size >= args.window:
      windows_count += 1
      #track progress
      if windows_count % 10000 == 0:
        print windows_count, "windows processed"
      # filter for missing data
      if all(i <= allowedN for i in genN):
        Midpos = int(words[1])- int(windows_size/2)
        CHR = words[0]
        outputNames.write('%s\t%s\n' % (CHR , Midpos))
        output = open(args.output+str(windows_count)+'.fasta', 'w') 
        for k, val in gend.iteritems():
          valP = ''.join(str(e) for e in val)
          output.write('>%s\n%s\n' % (k, valP))
      # reset the window:
      genotype = []
      gend = defaultdict(list)
      windows_size = 0
      genN = [0]*len(sampleNames)
    

    
    
      
# write last window if it meets the criteria
if windows_size >= args.window:
  windows_count += 1
  if all(i < allowedN for i in genN):
    output = open(args.output+str(windows_count)+'.fasta', 'w')
    for k, val in gend.iteritems():
      valP = ''.join(str(e) for e in val)
      output.write('>%s\n%s\n' % (k, valP))

datafile.close()
output.close()
