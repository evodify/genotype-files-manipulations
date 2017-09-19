#!/usr/bin/env python
'''
This script removes the masked sites from a tab file. The masked sites are provided in a BED file.

# test.tab:

Chr Pos Sample1 Saple2
scaffold_1  1 A C
scaffold_1  27 A A
scaffold_1  106 T A
scaffold_1  159 T A
scaffold_1  245 G G
scaffold_2  1 C C
scaffold_2  27 C C
scaffold_2  106 T C
scaffold_2  159 T C
scaffold_2  245 C C

# masked.bed:

CHR   start   end
scaffold_1      50   200
scaffold_2      115   159

# output.tab:

Chr Pos Sample1 Sample2
scaffold_1  1 A C
scaffold_1  27 A A
scaffold_1  245 G G
scaffold_2  1 C C
scaffold_2  27 C C
scaffold_2  106 T C
scaffold_2  245 C C

# command:

python mask_tab.py -i test.tab -m masked.bed -o output.tab

contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input BED file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-m', '--masked_intervals', help = 'file containing list of coordinates of masked regions in the BED format', type=str, required=True)
args = parser.parse_args()

############################# program #############################

maskFile = open(args.masked_intervals, "r")
header_scaf = maskFile.readline()
maskWords = maskFile.readline().split()
maskScaf = int(maskWords[0].split('_')[1])
maskStart = int(maskWords[1])
maskEnd = int(maskWords[2])

output = open(args.output, 'w')

counter = 0

with open(args.input) as datafile:
  header = datafile.readline()
  output.write("%s" % header)

  for line in datafile:
    words = line.split()
    inScaff = int(words[0].split('_')[1])
    inPos = int(words[1])

    # read masked coordinates until the position is in the window
    while  inScaff > maskScaf or (inScaff == maskScaf and inPos > maskEnd):
      maskWords = maskFile.readline().split()
      if maskWords == []:
        break
      else:
        maskScaf = int(maskWords[0].split('_')[1])
        maskStart = int(maskWords[1])
        maskEnd = int(maskWords[2])

    # skip if masked
    if inScaff == maskScaf and maskStart <= inPos and inPos <= maskEnd:
      continue
    else:
      output.write(line)

    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
output.close()
print('Done!')
