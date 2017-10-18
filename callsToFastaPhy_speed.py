#!/usr/bin/python2

"""
This script converts genotype calls file to FASTA and PHYLIP fast but consumes a lot of RAM. To convert calls to fasta/phy with low RAM consumption use callsToFastaPhy_RAM.py

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

$ python2 callsToFastaPhy_speed.py -i input.tab -o output -s "sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8"

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-f', '--fasta', help = 'name of the fasta output file', type=str, required=False)
parser.add_argument('-p', '--phylip', help = 'name of the phylip output file', type=str, required=False)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process (optional)', type=str, required=False)
args = parser.parse_args()

# check if any option was specified:
if not (args.fasta or args.phylip):
  raise IOError('Either -f or -p options need to be specified.')

# check if samples names are given and if all sample names are present in a header
sampleNames = calls.checkSampleNames(args.samples, args.input)

############################# program #############################

callsDF = calls.callsParser(args.input, sampleNames)

if args.fasta and args.phylip:
  outputFasta = open(args.fasta, 'w')
  outputPhy = open(args.phylip, 'w')

  NumberPos = len(callsDF.positions)
  NumberSamp = len(sampleNames)
  outputPhy.write(' %s %s\n' % (NumberSamp, NumberPos)) # make .phy header

  # write sample name into file
  for s in callsDF.names:
    outputFasta.write(">%s\n" % s)
    outputPhy.write("%s  " % s)
    
    seqCh = calls.chunks(callsDF[s], 100) # split sequence in multi-line fasta
    for c in seqCh:  # write sequence chunks
      cP = ''.join(str(i) for i in c)
      outputFasta.write("%s\n" % cP)
      outputPhy.write(cP)
    outputPhy.write("\n")
  
  outputFasta.close()
  outputPhy.close()
  
elif args.fasta:
  outputFasta = open(args.fasta, 'w')
  for s in callsDF.names:
    outputFasta.write(">%s\n" % s) # write sample name into file
    
    seqCh = calls.chunks(callsDF[s], 100) # split sequence in multi-line fasta
    for c in seqCh: # write sequence chunks
      cP = ''.join(str(i) for i in c)
      outputFasta.write("%s\n" % cP)
      
  outputFasta.close()
  
elif args.phylip:
  outputPhy = open(args.phylip, 'w')
  NumberPos = len(callsDF.positions)
  NumberSamp = len(sampleNames)
  outputPhy.write(' %s %s\n' % (NumberSamp, NumberPos))  # make .phy header
  
  for s in callsDF.names:
    outputPhy.write("%s  " % s)
    seqP = ''.join(str(i) for i in callsDF[s])
    outputPhy.write("%s\n" % seqP)

  outputPhy.close()
