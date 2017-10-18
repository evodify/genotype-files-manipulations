#! /usr/bin/python
'''
This script converts genotype calls file to ped and map files suitable for PLINK (http://zzz.bwh.harvard.edu/plink/).

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
chr_3   5   G   -   N   N   G   G   G   C   G
chr_4   1   G   -   N   N   G   G   G   C   G
chr_4   2   G   -   N   N   G   G   G   C   G

Diploid genotypes can be specified with /, e.g. A/T, T/C.

# output.ped
family1 sample1 0 0 0 0 A T C T 0 0 T T 0 0 A A C C 0 0 C C T T G G G C 0 0 T T - - - - - -
family1 sample2 0 0 0 0 0 0 T C - - T T C C A A C C 0 0 T T C C 0 0 C C 0 0 T T 0 0 0 0 0 0
family2 sample5 0 0 0 0 0 0 C C C C T T C C A A C C 0 0 C C T T 0 0 C C 0 0 T T G G G G G G
family2 sample6 0 0 0 0 0 0 0 0 C C T T C C A A C C 0 0 C C C T 0 0 0 0 0 0 T T G G G G G G

# output.map

1 1_1 0 1
1 1_2 0 2
1 1_3 0 3
1 1_4 0 4
1 1_6 0 6
2 2_1 0 1
2 2_2 0 2
2 2_3 0 3
2 2_4 0 4
2 2_5 0 5
3 3_1 0 1
3 3_2 0 2
3 3_3 0 3
3 3_4 0 4
3 3_5 0 5
4 4_1 0 1
4 4_2 0 2

# command

$ python calls_to_ped_map.py -i input.tab -o output -f "family1[sample1,sample2];family2[sample5,sample6]"

# An example command to convert to binary bed with PLINK:

$ plink --file output --make-bed --geno 0.999

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module
import re # to split input
import random # for randomization

############################# functions #############################
def pseudoPhasePED(gt):
  ''' Randomly splits heterozygouts '''
  # heterozygouts:
  ambR = ['A G', 'G A']
  ambY = ['T C', 'C T']
  ambM = ['A C', 'C A']
  ambK = ['G T', 'T G']
  ambS = ['G C', 'C G']
  ambW = ['A T', 'T A']
  delA = ['A -', '- A']
  delT = ['T -', '- T']
  delG = ['G -', '- G']
  delC = ['C -', '- C']
  phasedAlles = []
  for i in gt:
    if i == 'N':
      i = '0 0'
    elif i == 'A':
      i = 'A A'
    elif i == 'G':
      i = 'G G'
    elif i == 'C':
      i = 'C C'
    elif i == 'T':
      i = 'T T'
    elif i == '-' or i == '*' or i == "*/*":
      i = '- -'
    elif i == 'R' or i == 'A/G' or i == 'G/A':
      i = random.choice(ambR)
    elif i == 'Y' or i == 'T/C' or i == 'C/T':
      i = random.choice(ambY)
    elif i == 'M' or i == 'A/C' or i == 'C/A':
      i = random.choice(ambM)
    elif i == 'K' or i == 'G/T' or i == 'T/G':
      i = random.choice(ambK)
    elif i == 'S' or i == 'G/C' or i == 'C/G':
      i = random.choice(ambS)
    elif i == 'W' or i == 'A/T' or i == 'T/A':
      i = random.choice(ambW)
    elif i == 'A/*' or i == '*/A':
      i = random.choice(delA)
    elif i == 'T/*' or i == '*/T':
      i = rTndom.choice(delT)
    elif i == 'G/*' or i == '*/G':
      i = rGndom.choice(delG)
    elif i == 'C/*' or i == '*/C':
      i = rCndom.choice(delC)
    else:
      i = '0 0'
    phasedAlles.append(i)
  return phasedAlles

def familySampleCheck(family_samples, hap_dip_samples):
  ''' To check if samples specified in the options -f, -1n,-2n are the same'''
  if not set(family_samples) == set(hap_dip_samples):
    raise IOError('Sample names in -f differ from those specified in -1n, -2n')

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-f', '--family', help = 'Specify the family list in the format "family1[sample1,sample2];family2[sample5,sample6]"', type=str, required=True)
args = parser.parse_args()

# check Family list
familyNames = args.family
Fsamples = []
famDict = {}
for i in familyNames.strip("\"").split(";"):
  famName = i.split("[")[0]
  famSample = re.split("\[|\]|", i)[1]
  Fsamples.append(famSample.split(","))
  famDict[famName] = calls.checkSampleNames(famSample,args.input)
Fsamples = [i for sl in Fsamples for i in sl] # flat list

############################# program #############################

callsDF = calls.callsParser(args.input, Fsamples)

outputPED = open(args.output+'.ped', 'w')
outputMAP = open(args.output+'.map', 'w')

for i in range(len(callsDF[i])):
  # make map file
  snpsID = str(callsDF.chrmosomes[i]) + "_" + str(callsDF.positions[i])
  outputMAP.write("%s %s 0 %s\n" % (callsDF.chrmosomes[i], snpsID, callsDF.positions[i]))
outputMAP.close()

FamilyList = []
for i in range(len(callsDF.names)):
  ##create a family list
  for key, value in famDict.iteritems():
    if callsDF.names[i] in value:
       FamilyList.append(key)
      #print callsDF.names[i], f, vars()[f + "samples"]

  # make ped file
  seq = pseudoPhasePED(callsDF.sequences[i])
  seqP = ' '.join(str(s) for s in seq)

  outputPED.write("%s %s 0 0 0 0 %s\n" % (FamilyList[i], callsDF.names[i], seqP))
outputPED.close()

