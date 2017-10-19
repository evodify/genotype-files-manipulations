#!/usr/bin/python2

"""
This script removes sites with more than two alleles. Mononorphic sites are keept. Use removeMonomorphic.py to remove mononorphic sites.

# input:

CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   1   A   W   N   N   A   N   N   N   N
chr_1   2   C   Y   Y   N   C   C   N   C   N
chr_1   3   C   N   -   N   C   C   C   C   C
chr_1   4   C   A   T   N   T   T   T   T   T
chr_1   6   G   Y   C   N   C   C   C   C   C
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
chr_2   1   A   A   A   N   A   A   A   A   A
chr_2   2   C   C   C   N   C   C   C   C   C
chr_2   3   C   N   N   N   N   N   N   N   N
chr_2   4   C   C   T   C   C   C   C   C   C
chr_2   5   T   T   C   T   Y   T   Y   T   T
chr_3   1   G   G   N   N   G   N   N   N   N
chr_3   2   C   S   C   N   C   C   N   C   N
chr_3   4   N   T   T   N   T   T   T   T   N

# command:

$ python keep_biallelic.py -i inputfile -o outputfile

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

with open(args.input) as datafile:

  print('Creating the output file...')
  header_line = datafile.readline()
  fileoutput = open(args.output, 'w')
  fileoutput.write(header_line)

  for line in datafile:
    words = line.split()
    chr_pos = words[0:2]
    genotypes = words[2:]

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

    # filter by alleles number
    if set(['N']).issuperset(genotypes):  # skip all Ns sites
      continue
    elif set(['A', 'C', 'M', 'N']).issuperset(genotypes) or set(['A', 'T', 'W', 'N']).issuperset(genotypes) or set(['A', 'G', 'R', 'N']).issuperset(genotypes) or set(['T', 'C', 'Y', 'N']).issuperset(genotypes) or set(['T', 'G', 'K', 'N']).issuperset(genotypes) or set(['G', 'C', 'S', 'N']).issuperset(genotypes):
      fileoutput.write(line)

fileoutput.close()
datafile.close()
print('Done!')
