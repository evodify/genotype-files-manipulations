#!/usr/bin/python2

"""
This script phases the sequences by random split of heterozygous sites.

# input:

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

# output:

CHROM   POS REF_1   REF_2   sample1_1   sample1_2   sample2_1   sample2_2   sample3_1   sample3_2   sample4_1   sample4_2   sample5_1   sample5_2   sample6_1   sample6_2   sample7_1   sample7_2   sample8_1   sample8_2 
chr_1   1   A   A   A   T   N   N   N   N   A   A   N   N   N   N   N   N   N   N
chr_1   2   C   C   C   T   T   C   N   N   C   C   C   C   N   N   C   C   N   N
chr_1   3   C   C   N   N   -   -   N   N   C   C   C   C   C   C   C   C   C   C
chr_1   4   T   T   T   T   T   T   N   N   T   T   T   T   T   T   T   T   T   T
chr_1   6   C   C   N   N   C   C   N   N   C   C   C   C   C   C   C   C   C   C
chr_2   1   A   A   A   A   A   A   N   N   A   A   A   A   A   A   A   A   A   A
chr_2   2   C   C   C   C   C   C   N   N   C   C   C   C   C   C   C   C   C   C
chr_2   3   C   C   N   N   N   N   N   N   N   N   N   N   N   N   N   N   N   N
chr_2   4   C   C   C   C   T   T   C   C   C   C   C   C   C   C   C   C   C   C
chr_2   5   T   T   T   T   C   C   T   T   T   C   T   T   C   T   T   T   T   T
chr_3   1   G   G   G   G   N   N   N   N   G   G   N   N   N   N   N   N   N   N
chr_3   2   C   C   C   G   C   C   N   N   C   C   C   C   N   N   C   C   N   N
chr_3   3   N   N   N   N   N   N   N   N   N   N   N   N   N   N   N   N   N   N
chr_3   4   N   N   T   T   T   T   N   N   T   T   T   T   T   T   T   T   N   N
chr_3   5   G   G   -   -   N   N   N   N   G   G   G   G   G   G   C   C   G   G
chr_4   1   G   G   -   -   N   N   N   N   G   G   G   G   G   G   C   C   G   G
chr_4   2   G   G   -   -   N   N   N   N   G   G   G   G   G   G   C   C   G   G

# command:

$ python pseudoPhasingHetero.py -i inputfile -o outputfile

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
args = parser.parse_args()

############################# program #############################


counter = 0

print('Opening the file...\n')

with open(args.input) as datafile:
  header_words = datafile.readline().split()

  attrib = header_words[0:2]
  colbeg = '\t'.join(str(e) for e in attrib)

  # create slitted column names
  ids = header_words[2:]
  idsheader = []
  for i in ids:
    idsheader.append("%s_1\t%s_2 " % (i, i))
  header = '\t'.join(str(e) for e in idsheader)

  print('Creating the output file...\n')
  outname = '%s' % args.output
  outfile = open(outname, 'w')
  outfile.write('%s\t%s\n' % (colbeg, header))  

  print('Randomizing alleles in heterozygous genotypes ...\n')

  for line in datafile:
    words = line.split()
    chr_pos = words[0:2]
    alleles = words[2:]

    # check if all genotypes are correct
    calls.if_all_gt_correct(alleles, line)

    # phase
    assignAlles = calls.pseufoPhase(alleles)
 
    chromPosP = '\t'.join(str(e) for e in chr_pos)
    assignAllesP = '\t'.join(str(e) for e in assignAlles)
    outfile.write("%s\t%s\n" % (chromPosP, assignAllesP))
    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
outfile.close()
print('\nDone!')