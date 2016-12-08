#!/usr/bin/python2

"""
This script converts genotype calls file to FASTA and PHYLIP.

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

# output.fasta:

>sample1
WYNTNACNCTGSNT---
>sample2
NY-TCACNTCNCNTNNN
>sample3
NNNNNNNNCTNNNNNNN
>sample4
ACCTCACNCYGCNTGGG
>sample5
NCCTCACNCTNCNTGGG
>sample6
NNCTCACNCYNNNTGGG
>sample7
NCCTCACNCTNCNTCCC
>sample8
NNCTCACNCTNNNNGGG

# output.phy:

 8 17
sample1 WYNTNACNCTGSNT---
sample2 NY-TCACNTCNCNTNNN
sample3 NNNNNNNNCTNNNNNNN
sample4 ACCTCACNCYGCNTGGG
sample5 NCCTCACNCTNCNTGGG
sample6 NNCTCACNCYNNNTGGG
sample7 NCCTCACNCTNCNTCCC
sample8 NNCTCACNCTNNNNGGG

# command:

$ python2 callsToFastaPhy.py -i input.tab -o output -s "sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8"

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process (optional)', type=str, required=False)
args = parser.parse_args()

# check if samples names are given and if all sample names are present in a header
sampleNames = calls.checkSampleNames(args.samples, args.input)

############################# program #############################

# count number of samples and positions for .phy header
NumberPos = calls.countPositions(args.input)
NumberSamp = len(sampleNames)

outputFasta = open(args.output+'.fasta', 'w')
outputPhy = open(args.output+'.phy', 'w')

# make .phy header 
outputPhy.write(' %s %s\n' % (NumberSamp, NumberPos))

# process one sample per time to reduce RAM usage
for sample in sampleNames:

  # write sample name into file
  outputFasta.write(">%s\n" % sample)
  outputPhy.write("%s  " % sample)

  fastaLim = 0 # counter to split sequence in multi-line fasta

  with open(args.input) as datafile:
    header_words = datafile.readline().split()

    # index a sample
    sampCol = calls.indexSamples([sample], header_words)

    for line in datafile:
      words = line.split()

      genotype = calls.selectSamples(sampCol, words)

      # output only single nucleotide genotypes, insertions are replaced with N.
      if len(genotype) == 1:
        outputFasta.write(genotype[0])
        outputPhy.write(genotype[0])
      else:
        outputFasta.write('N')
        outputPhy.write('N')

      # to split sequence in multi-line fasta
      fastaLim += 1
      if fastaLim == 100:
        outputFasta.write("\n")
        fastaLim = 0

  outputFasta.write("\n")
  outputPhy.write("\n")
  print sample, 'processed'

datafile.close()
outputFasta.close()
outputPhy.close()
