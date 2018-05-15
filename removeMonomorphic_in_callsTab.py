#!/usr/bin/python2

"""
This script removes monomorphic positions, i.e. keeps only SNPs.

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

CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   1   A   W   N   N   A   N   N   N   N
chr_1   2   C   Y   Y   N   C   C   N   C   N
chr_1   3   C   N   -   N   C   C   C   C   C
chr_2   4   C   C   T   C   C   C   C   C   C
chr_2   5   T   T   C   T   Y   T   Y   T   T
chr_3   2   C   S   C   N   C   C   N   C   N
chr_3   5   G   -   N   N   G   G   G   C   G
chr_4   1   G   -   N   N   G   G   G   C   G
chr_4   2   G   -   N   N   G   G   G   C   G

# command:

$ python removeMonomorphic.py -i inputfile -o outputfile

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""
############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
args = parser.parse_args()

############################# program #############################

counter = 0

print('Opening the file...')

fileoutput = open(args.output, 'w')

with open(args.input) as datafile:
  header_line = datafile.readline()
  fileoutput.write(header_line)

  print('Processing the file ...')

  for line in datafile:
    words = line.split()
    GT = words[2:]

    if '/' in GT[1]:
      GT = calls.twoToOne(GT)

    # output the line if it is polymorphic
    if calls.is_polymorphic(GT):
      fileoutput.write(line)

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
fileoutput.close()

print('Done!')
