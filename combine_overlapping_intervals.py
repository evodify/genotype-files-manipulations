#!/usr/bin/env python
'''
This script combines overlapping genetic intervals in the bed format.


input.file:

CHR	start	end
scaffold_1	1704	2714
scaffold_1	2704	3704
scaffold_1	3804	4704
scaffold_1	4904	5704
scaffold_1	6704	7704
scaffold_1	7704	8704
scaffold_1	8704	9704
scaffold_1	9704	10704

input.file:


CHR	start	end
scaffold_1	1604	5804
scaffold_1	6704	10704

command:

python combine_overlapping_intervals.py -i input.file -o output.file -v 100


contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-v', '--overhang', help = 'overhang size', type=int, required=False)

args = parser.parse_args()

############################# program #############################

if not args.overhang:
  overhang = 0
else: 
  overhang = args.overhang

output = open(args.output, 'w')

counter = 0

startPrevious = 0
endPrevious = 0

with open(args.input) as datafile:
  header = datafile.readline()
  output.write("%s" % header)

  for line in datafile:
    words = line.split()
    chrom = words[0]
    startNext = int(words[1])
    endNext = int(words[2])

    if (startNext-overhang) <= (endPrevious+overhang):
      endPrevious = endNext
    elif startPrevious == 0:
      startPrevious = startNext
      endPrevious = endNext
    else:
      output.write("%s\t%s\t%s\n" % (chrom, startPrevious-overhang, endPrevious+overhang))
      startPrevious = startNext
      endPrevious = endNext

    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"
      
  output.write("%s\t%s\t%s\n" % (chrom, startPrevious, endPrevious))
  
datafile.close()
output.close()
print('Done!')
