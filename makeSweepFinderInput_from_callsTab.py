#!/usr/bin/python2

"""
This script makes an input file for [SweepFinder](http://people.binf.ku.dk/rasmus/webpage/sf.html). To make such file information on the outgroup/ancestral sequence is required.

Note! Chromosome number in the first column must be separated by _.
For example, chr_1 - correct, chr1 - incorrect.

# input:

#CHROM  POS sample1 sample2 sample3 sample4
chr_1  113 G   G   G   N
chr_1  117 C   C   N   C
chr_1  122 C   C   N   C
chr_2  137 A   A   T   N
chr_2  139 T   T   T   N
chr_2  148 A   A   T   N
chr_13  161 C   T   C   N
chr_13  170 C   T   C   N
chr_13  174 A   A   A   N
chr_X  104 A   A   A   N
chr_X  204 G   G   G   N
chr_X  574 A   A   A   N
chr_Un1  174 A   A   A   N
chr_Un2  104 A   A   A   N
chr_Un4  204 G   G   G   N
chr_Un7  574 A   A   A   N

# ancestral:

#CHROM  POS Co_ancestor
chr_1  113 G
chr_1  117 N
chr_1  122 C
chr_2  137 A
chr_2  139 T
chr_2  148 T
chr_13  161 S
chr_13  170 C
chr_13  174 A
chr_X  104 A
chr_X  204 A
chr_X  574 A
chr_Un1  174 T
chr_Un2  104 G
chr_Un4  204 C
chr_Un7  574 C

# fai:

chr_1      196        12      80      81
chr_2      1410        198        80      81
chr_13      1500        341        80      81
chr_X      10000        301        80      81
chr_Un1  20000        101        80      81
chr_Un2  30000        301        80      81
chr_Un4  40000        401        80      81
chr_Un7  50000        901        80      81

# output (There are also separate output files for each chromosome.):

position	x	n	folded
117	3	3	1
333	1	3	0
344	2	3	0
1767	1	3	0
1776	1	3	0
3310	3	3	0


# command:

$ python makeSweepFinderInput_from_callsTab.py -i test.tab -o test -f test.fai -a ancestral.test -N 1 -m chr_X

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""
############################# modules #############################

import calls  # my custom module
import collections
import re
import random

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help='name of the input file', type=str, required=True)
parser.add_argument('-a', '--ancestor', help='name of the outgroup/ancestral sequence file', type=str, required=True)
parser.add_argument('-f', '--fai', help='name of the fasta.fai file', type=str, required=True)
parser.add_argument('-o', '--output', help='name of the output file', type=str, required=True)
parser.add_argument('-N', '--missing', help='number of allowed Ns', type=int, required=True)
parser.add_argument('-m', '--major_chromosomes', help='name of the last major chromosome', type=str, required=True)
parser.add_argument('-s', '--samples', help='column names of the samples to process (optional)', type=str,
                    required=False)
args = parser.parse_args()

# check if samples names are given and if all sample names are present in a header
sampleNames = calls.checkSampleNames(args.samples, args.input)

# check missing data settings
Ns = calls.checkMissing(args.missing)

# check the name of the last major chromosome
mChr = args.major_chromosomes.split('_')[1]
try: # for sex chromosomes
    mChr = int(mChr)
except Exception:
    mChr = ord(mChr)

############################# program #############################

counter = 0

noBreak = bool(True)
Ns = args.missing
fai = open(args.fai, 'r')
fai_words = fai.readline().split()
fai_ch = fai_words[0]
fai_start1 = 0
fai_start2 = int(fai_words[1])

output = open(args.output, 'w')
output.write("position\tx\tn\tfolded\n")
outputChr = open(fai_ch + '_' + args.output, 'w')
outputChr.write("position\tx\tn\tfolded\n")

ancestFile = open(args.ancestor, 'r')
ancest_header = ancestFile.readline()
words2 = ancestFile.readline().split()
ancest_chr_pos = words2[0:2]
ancest_ch = ancest_chr_pos[0].split('_')[1]
try: # for sex chromosomes
  ancest_ch = int(ancest_ch)
except Exception:
  ancest_ch = ord(ancest_ch)


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
        ch = words[0].split('_')[1]
        try: # for sex chromosomes
          ch = int(ch)
        except Exception:
            try:
                ch = ord(ch)
            except Exception:
                ch = ch
        pos = int(words[1])

        # track progress
        counter += 1
        if counter % 1000000 == 0:
            print str(counter), "lines processed"

        # select alleles
        alleles = calls.selectSamples(sampCol, words)

        # check if one- or two-character code
        if any(["/" in gt for gt in alleles]):
            alleles = calls.twoToOne(alleles)

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
                ancest_ch = ancest_chr_pos[0].split('_')[1]
                try: # for sex chromosomes
                  ancest_ch = int(ancest_ch)
                except Exception:
                    try:
                        ancest_ch = ord(ancest_ch)
                    except Exception:
                        ancest_ch = ancest_ch
                ancest_pos = int(ancest_chr_pos[1])
                ancest = words2[2]
                if ancest in 'RYMKSW':
                    ancest = calls.OneToTwo(ancest)[0].split('/')

        # find overlap with fai file to define chromosome borders
        if ch <= mChr:  # major chromosomes that will be split
            while chr_pos[0] != fai_ch:
                fai_words = fai.readline().split()
                if fai_words == []:
                    break
                else:
                    fai_ch = fai_words[0]
                    fai_start1 = fai_start1 + fai_start2
                    fai_start2 = int(fai_words[1])
                    outputChr = open(fai_ch + '_' + args.output, 'w')
                    outputChr.write("position\tx\tn\tfolded\n")
            posP = pos + fai_start1
        #  break after major scaffolds
        else:
            break

        # polarize alleles
        n = len(AllallesNoN)
        al1 = numAl[0][0]
        x1 = numAl[0][1]
        f = 0
        if len(numAl) > 2:  # skip non-biallelic
            continue
        elif chr_pos == ancest_chr_pos:  # folded
            if len(numAl) == 1:  # fixed
                if ancest == 'N':
                    f = 1
                    x = random.choice([x1, 0])
                elif al1 in ancest:
                    x = 0
                else:
                    x = n
            elif len(numAl) == 2:  # biallelic
                al2 = numAl[1][0]
                x2 = numAl[1][1]
                if (al1 in ancest) and (al2 not in ancest):
                    x = x2
                elif (al2 in ancest) and (al1 not in ancest):
                    #print ancest, al2, al1
                    x = x1
                else:

                    f = 1
                    x = random.choice([x1, x2])
        else:  # unfolded
            if len(numAl) == 1:  # fixed
                f = 1
                x = random.choice([x1, 0])
            elif len(numAl) == 2:
                x2 = numAl[1][1]
                f = 1
                x = random.choice([x1, x2])
        #print ancest, posP, numAl, x, n, f, line

        if x == 0 and f == 0: # skip sites with fixed ancestral alleles
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
