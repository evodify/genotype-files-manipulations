#!/usr/bin/env python
'''
This script converts a tab-delimited file to a bed file.

input.file:

CHR	pos
scaffold_1      1
scaffold_1      2
scaffold_1      3
scaffold_1      4
scaffold_1      10
scaffold_1      16
scaffold_1      17
scaffold_2      1
scaffold_2      2
scaffold_2      10

output.file:

Chr Start   End
scaffold_1  1   17
scaffold_2  1   2
scaffold_2  10  10

command:

python2 tabToBED.py -i test.tab -o test.bed -c 5

contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-c', '--collapse', help = 'size between intervals to collapse ', type=int, required=False)

args = parser.parse_args()

############################# program #############################

if not args.collapse:
  collapse = 0
else:
  collapse = args.collapse

output = open(args.output, 'w')

counter = 0
chromP = []
posPrevious = 0
endPrevious = 0

with open(args.input) as datafile:
  header = datafile.readline()
  output.write("Chr\tStart\tEnd\n")
  for line in datafile:
    words = line.split()
    chrom = words[0]
    pos = int(words[1])

    if chromP == []:
      posStart = pos
      chromP = chrom
      posPrevious = pos
      chromPrevious = chrom
    elif chrom == chromPrevious and (posPrevious+collapse) >= (pos-1):
      posPrevious = pos
      chromPrevious = chrom
    else:
      output.write("%s\t%s\t%s\n" % (chromP, posStart, posPrevious))
      posStart = pos
      chromP = chrom
      posPrevious = pos
      chromPrevious = chrom

    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

  output.write("%s\t%s\t%s\n" % (chromP, posStart, posPrevious))

datafile.close()
output.close()
print('Done!')
