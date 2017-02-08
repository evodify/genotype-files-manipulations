#!/usr/bin/env python

'''
This script calculates number of positions with missing data (Ns) using sliding window approach.

Example input:

CHROM   POS SomeValues
chr_1   1   2456
chr_1   2   36
chr_1   3   346
chr_1   4   36
chr_2   1   36
chr_2   2   2461
chr_2   3   6
chr_2   4   70
chr_2   5   2464
chr_3   1   46
chr_3   2   2466
chr_3   3   36
chr_3   4   6
chr_3   5   6

Example fasta.fai:

chr_1      19624517
chr_2      14106692
chr_3      15060676

Example output:

position    SomeValues
1   2456
2   36
3   346
4   36
19624518    36
19624519    2461
19624520    6
19624521    70
19624522    2464
33731210    46
33731211    2466
33731212    36
33731213    6
33731214    6

command:

$ python mergeChrPos.py -i input.tab -F fasta.fai -o output.file

contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''
############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-f', '--fai', help = 'name of the fasta.fai file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
args = parser.parse_args()

############################# program #############################

counter = 0

print('Opening files...')

fai = open(args.fai, 'r')
fai_words = fai.readline().split()
fai_ch = fai_words[0]
fai_start1 = 0
fai_start2 = int(fai_words[1])

output = open(args.output, 'w')

print('Opening the file...')
with open(args.input) as datafile:
  header_line = datafile.readline().split()
  header = header_line[2:]
  headerP = '\t'.join(str(e) for e in header)
  output.write("position\t%s\n" % headerP)
  for line in datafile:
    words = line.split()
    chr = words[0]
    pos = float(words[1])
    values = words[2:]
    valuesP = '\t'.join(str(e) for e in values)

    while chr != fai_ch:
      fai_words = fai.readline().split()
      if not fai_words:  # end of fai file
        raise Exception("%s is not found in %s.\nExecution stopped!" % (chr, args.fai))
      else:
        fai_ch = fai_words[0]
        fai_start1 = fai_start1+fai_start2
        fai_start2 = float(fai_words[1])
    posP = pos+fai_start1

    output.write("%s\t%s\n" % (posP, valuesP))

    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
output.close()
print('Done!')
 
