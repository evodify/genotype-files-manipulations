#!/usr/bin/python2

"""
This script removes insertions of longer than 1 bp and replaces deletions of 1 bp marked as "*" with "-".

# input:

CHROM   POS REF sample1 sample2 sample3
chr_1   1   A   G   N 
chr_1   2   CCATTAGAT   A   C
chr_1   3   C   N   *

# output:

CHROM   POS REF sample1 sample2 sample3
chr_1   1   A   G   N
chr_1   2   N   A   C
chr_1   3   C   N   -


# command:

$ python remove_Isertions.py -i inputfile -o outputfile

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

print('Opening the file...')

fileoutput = open(args.output, 'w')

with open(args.input) as datafile:
  header_line = datafile.readline()
  fileoutput.write(header_line)

  print('Processing the file ...')

  for line in datafile:
    words = line.split()
    ChrPos = words[0:2]
    GT = words[2:]

    # output the line if it is polymorphic
    for i in range(0,len(GT)):
      if len(GT[i])!= 1:
        GT[i] = "N"
      elif GT[i]=="*":
        GT[i] = "-"
    ChrPosP = '\t'.join(str(e) for e in ChrPos)
    GTP = '\t'.join(str(e) for e in GT)

    fileoutput.write("%s\t%s\n" % (ChrPosP, GTP))

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
fileoutput.close()

print('Done!')
