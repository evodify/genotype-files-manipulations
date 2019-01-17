#! /usr/bin/python
'''
This script converts genotype calls file to ped and map files 
suitable for PLINK (http://zzz.bwh.harvard.edu/plink/).

Note! If your file is large, this script may use a lot of RAM.

# input.tab
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
chr_2   5   T   T   C   T   C   T   Y   T   T
chr_3   1   G   G   N   N   G   N   N   N   N
chr_3   2   C   S   C   N   C   C   N   C   N
chr_3   3   N   N   N   N   N   N   N   N   N
chr_3   4   N   T   T   N   T   T   T   T   N

Diploid genotypes can be specified with /, e.g. A/T, T/C.

If the data is phased, genotypes need to be separated by |, e.g C|C	A|C

# output.ped
family1 sample1 0 0 0 0 A T C T 0 0 T T 0 0 A A C C 0 0 C C T T G G G C 0 0 T T 
family1 sample2 0 0 0 0 0 0 T C - - T T C C A A C C 0 0 T T C C 0 0 C C 0 0 T T
family2 sample5 0 0 0 0 0 0 C C C C T T C C A A C C 0 0 C C T T 0 0 C C 0 0 T T
family2 sample6 0 0 0 0 0 0 0 0 C C T T C C A A C C 0 0 C C C T 0 0 0 0 0 0 T T

# output.map

chr_1 chr_1_1 0 1
chr_1 chr_1_2 0 2
chr_1 chr_1_3 0 3
chr_1 chr_1_4 0 4
chr_1 chr_1_6 0 6
chr_2 chr_2_1 0 1
chr_2 chr_2_2 0 2
chr_2 chr_2_3 0 3
chr_2 chr_2_4 0 4
chr_2 chr_2_5 0 5
chr_3 chr_3_1 0 1
chr_3 chr_3_2 0 2
chr_3 chr_3_3 0 3
chr_3 chr_3_4 0 4

# command

$ python calls_to_ped_map.py \
-i input.tab \
-o output \
-f "family1[sample1,sample2];family2[sample5,sample6]"

# An example command to convert to binary bed with PLINK:

$ plink --file output --make-bed --geno 0.999

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module
import re # to split input
import random # for randomization

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument(
  '-i',
  '--input',
  help = 'name of the input file',
  type = str,
  required = True)
parser.add_argument(
  '-o',
  '--output',
  help = 'name of the output file',
  type = str,
  required = True)
parser.add_argument(
  '-f',
  '--family',
  help = 'Specify the family list in the format \
          "family1[sample1,sample2];family2[sample5,sample6]"',
  type=str,
  required=True)
args = parser.parse_args()

# check Family list
familyNames = args.family
Fsamples = []
famDict = {}
for i in familyNames.strip("\"").split(";"):
  famName = i.split("[")[0]
  famSample = re.split("\[|\]|", i)[1]
  Fsamples.append(famSample.split(","))
  famDict[famName] = calls.checkSampleNames(famSample, args.input)
Fsamples = calls.flattenList(Fsamples)

############################# program #############################

callsDF = calls.callsParser(args.input, Fsamples)

outputPED = open(args.output+'.ped', 'w')
outputMAP = open(args.output+'.map', 'w')

for i in range(len(callsDF.positions)):
  # make map file
  snpsID = str(callsDF.chrmosomes[i]) + "_" + str(callsDF.positions[i])
  cM = float(callsDF.positions[i])*0.000001
  outputMAP.write("%s %s %s %s\n" %
                 (callsDF.chrmosomes[i], snpsID, cM,callsDF.positions[i]))
outputMAP.close()

FamilyList = []
for i in range(len(callsDF.names)):
  ##create a family list
  for key, value in famDict.iteritems():
    if callsDF.names[i] in value:
       FamilyList.append(key)

  # make ped file
  if all("|" in gt for gt in callsDF.sequences[i]):
    seq = calls.phasePED(callsDF.sequences[i])
  else:
    seq = calls.pseudoPhasePED(callsDF.sequences[i])
  seqP = ' '.join(str(s) for s in seq)

  outputPED.write("%s %s 0 0 0 0 %s\n" % (FamilyList[i], callsDF.names[i], seqP))
outputPED.close()

