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
family1 sample1 0   0   0   0   A T T C N N T T N N A A C C N N C C T T G G C G N N T T - - - - - -
family1 sample2 0   0   0   0   N N T C - - T T C C A A C C N N T T C C N N C C N N T T N N N N N N
family2 sample5 0   0   0   0   N N C C C C T T C C A A C C N N C C T T N N C C N N T T G G G G G G
family2 sample6 0   0   0   0   N N N N C C T T C C A A C C N N C C Y Y N N N N N N T T G G G G G G

# output.map

1   chr_1_1 0   1
1   chr_1_2 0   2
1   chr_1_3 0   3
1   chr_1_4 0   4
1   chr_1_6 0   6
2   chr_2_1 0   1
2   chr_2_2 0   2
2   chr_2_3 0   3
2   chr_2_4 0   4
2   chr_2_5 0   5
3   chr_3_1 0   1
3   chr_3_2 0   2
3   chr_3_3 0   3
3   chr_3_4 0   4
3   chr_3_5 0   5
4   chr_4_1 0   1
4   chr_4_2 0   2

# command

$ python calls_to_ped_map.py -i input.tab -o output -f "family1[sample1,sample2];family2[sample5,sample6]" -2n "sample1,sample2" -1n "sample5,sample6"

# An example command to convert to binary bed with PLINK:

$ plink --file output --make-bed --geno 0.999

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module
import pandas as pd
import random # for randomization
import re # to split input

############################# functions ###########################

def pseudoPhaseSpace(gt):
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
  if gt == 'N':
    gt = 'N N'
  elif gt == 'A':
    gt = 'A A'
  elif gt == 'G':
    gt = 'G G'
  elif gt == 'C':
    gt = 'C C'
  elif gt == 'T':
    gt = 'T T'
  elif gt == '-' or gt == '*' or gt == "*/*":
    gt = '- -'
  elif gt == 'R' or gt == 'A/G' or gt == 'G/A':
    gt = random.choice(ambR)
  elif gt == 'Y' or gt == 'T/C' or gt == 'C/T':
    gt = random.choice(ambY)
  elif gt == 'M' or gt == 'A/C' or gt == 'C/A':
    gt = random.choice(ambM)
  elif gt == 'K' or gt == 'G/T' or gt == 'T/G':
    gt = random.choice(ambK)
  elif gt == 'S' or gt == 'G/C' or gt == 'C/G':
    gt = random.choice(ambS)
  elif gt == 'W' or gt == 'A/T' or gt == 'T/A':
    gt = random.choice(ambW)
  elif gt == 'A/*' or gt == '*/A':
    gt = random.choice(delA)
  elif gt == 'T/*' or gt == '*/T':
    gt = rTndom.choice(delT)
  elif gt == 'G/*' or gt == '*/G':
    gt = rGndom.choice(delG)
  elif gt == 'C/*' or gt == '*/C':
    gt = rCndom.choice(delC)
  else:
    gt = 'N N'
  return gt

def familySampleCheck(family_samples, hap_dip_samples):
  ''' To check if samples specified in the options -f, -1n,-2n are the same'''
  if not set(family_samples) == set(hap_dip_samples):
    raise IOError('Sample names in -f differ from those specified in -1n, -2n')

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-f', '--family', help = 'Specify the family list in the format "family1[sample1,sample2];family2[sample5,sample6]"', type=str, required=True)
parser.add_argument('-1n', '--haploid', help = 'Specify the ulations which are haploid in the format "sample1,sample2". Samples are considered diploid by default.', type=str, required=False)
parser.add_argument('-2n', '--diploid', help = 'If you specified haploid samples, also specify the ulations which are diploid in the format "sample3,sample4". Samples that are not specified will not be output.', type=str, required=False)
args = parser.parse_args()

# check Family list
familyNames = args.family
families = []
Fsamples = []
for i in familyNames.strip("\"").split(";"):
  famName = i.split("[")[0]
  famSample = re.split("\[|\]|", i)[1]
  families.append(famName)
  Fsamples.append(famSample.split(","))
  vars()[famName + "samples"] = calls.checkSampleNames(famSample,args.input)
Fsamples = [i for sl in Fsamples for i in sl] # flat list

# check sample names
if args.haploid and args.diploid:
  Hnames = calls.checkSampleNames(args.haploid, args.input)
  Dnames = calls.checkSampleNames(args.diploid, args.input)
  familySampleCheck(Fsamples, Dnames+Hnames)
elif args.haploid :
  Hnames = calls.checkSampleNames(args.haploid, args.input)
  familySampleCheck(Fsamples, Hnames)
elif args.diploid:
  Dnames = calls.checkSampleNames(args.diploid, args.input)
  familySampleCheck(Fsamples, Dnames)
else:
  Dnames = Fsamples
  print "No samples were specified.\nAll samples will be output and processed as diploids..."

############################# program #############################

outputName = args.output
callsDF = pd.read_csv(args.input, delim_whitespace=True)

#####################
# create a map file #
#####################

mapDF = pd.DataFrame({ '0Chr' : callsDF[callsDF.columns[0]].str.split('_',1).str[1],
                 '1snpsID' : callsDF[callsDF.columns[0]] + "_" + callsDF[callsDF.columns[1]].map(str),
                 '2dist' : 0,
                 '3Pos' : callsDF[callsDF.columns[1]]})
# output map
mapDF.to_csv(outputName+'.map', sep='\t', header=False, index=False)
del mapDF # clear mapDF from the memory

#####################
# create a ped file #
#####################

if args.haploid: # duplicate haploid genotypes
  for h in Hnames:
    callsDF.loc[:,(h)] = callsDF[h]+" "+callsDF[h]
    
### split diploid genotypes:
if args.diploid:
  for d in Dnames:
    for i in range(len(callsDF[d])):
      #print len(callsDF[d]), d, callsDF[d][i]
      if callsDF[d][i] in "ATGCN-":
        callsDF.loc[i,d] = callsDF[d][i]+" "+callsDF[d][i]
      else:
        callsDF.loc[i,d] = pseudoPhaseSpace(callsDF[d][i])
      #print len(callsDF[d]), d, callsDF[d][i]

pedDF = callsDF[Fsamples].T # transpose and drop CHR,POS

del callsDF # clear callsDF from the memory

#create a family list
FamilyList = []
for s in list(pedDF.index):
  for f in families:
    if s in vars()[f + "samples"]:
      FamilyList.append(f)
      #print s,f,vars()[f + "samples"]

# add the necessary columns
pedDF.insert(0, "Family", FamilyList) # insert Family ID column
pedDF.insert(1, "Sample", list(pedDF.index)) # insert Sample ID column
pedDF.insert(2, "PaternID", 0) # insert Paternal ID column
pedDF.insert(3, "MaternID", 0) # insert Maternal ID column
pedDF.insert(4, "Sex", 0) # insert Sex column
pedDF.insert(5, "Aff", 0) # insert Affection column

# output ped
pedDF.to_csv(outputName+'.ped', sep='\t', header=False, index=False)
