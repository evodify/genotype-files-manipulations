#! /usr/bin/env python

"""
This script merges two tab files by the overlapping positions.
All the positions of the input2 file that do not exist in the reference input will not be output. While position that exist in the reference but do not exist in the input2 will be represented by Ns.

# reference_input:

CHROM   POS sample1 sample2
chr_1   1   A   W
chr_1   2   C   Y
chr_1   3   C   N
chr_1   4   T   T
chr_1   6   C   C
chr_2   1   A   A
chr_2   2   C   C
chr_2   3   C   N
chr_2   4   C   C
chr_2   5   T   T
chr_3   1   G   G
chr_3   2   C   S
chr_3   3   N   N
chr_3   4   N   T
chr_3   5   G   -

# input2:

CHROM   POS sample11    sample21
chr_1   1   N   N
chr_1   2   C   N
chr_1   3   C   C
chr_1   4   T   T
chr_1   7   T   T
chr_2   1   A   A
chr_2   2   C   C
chr_2   3   N   N
chr_2   4   C   C
chr_2   5   T   T
chr_3   1   N   N
chr_3   2   C   N
chr_3   3   N   N
chr_3   4   T   N
chr_3   5   C   G

# output:

CHROM   POS sample1 sample2 sample11    sample21
chr_1   1   A   W   N   N
chr_1   2   C   Y   C   N
chr_1   3   C   N   C   C
chr_1   4   T   T   T   T
chr_1   6   C   C   N   N
chr_2   1   A   A   A   A
chr_2   2   C   C   C   C
chr_2   3   C   N   N   N
chr_2   4   C   C   C   C
chr_2   5   T   T   T   T
chr_3   1   G   G   N   N
chr_3   2   C   S   C   N
chr_3   3   N   N   N   N
chr_3   4   N   T   T   N
chr_3   5   G   -   C   G


# command:

$ python mergeTabFiles.py -s reference_input -l input2 -o output

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""
############################# modules #############################

import argparse, sys

############################# options #############################


class MyParser(argparse.ArgumentParser): 
   def error(self, message):
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

parser = MyParser()
parser.add_argument('-r', '--reference_input', help = 'name of the reference input file', type=str, required=True)
parser.add_argument('-i', '--large_input', help = 'name of the larger input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)

args = parser.parse_args()

############################# program #############################

counter = 0

# read the header of the larger file
linput = open(args.large_input, 'r')
linput_words = linput.readline().split()
linput_ch = linput_words[0]
linput_pos = linput_words[1]
linput_gt = linput_words[2:]

output = open(args.output, 'w')

print('Opening the file...')
with open(args.reference_input) as datafile:
  header_line = datafile.readline().split()
  print('Creating the output file...')
  headerP = '\t'.join(str(e) for e in header_line+linput_gt)
  output.write('%s\n' % headerP)

  # read the second line of the larger file
  linput_words = linput.readline().split()
  linput_ch = int(linput_words[0].split('_')[1])
  linput_pos = int(linput_words[1])
  linput_gt = linput_words[2:]

  for line in datafile:
    words = line.split()
    ch = int(words[0].split('_')[1])
    pos = int(words[1])

    # find overlap
    while (ch > linput_ch) or (ch == linput_ch and pos > linput_pos):
      linput_words = linput.readline().split()
      if linput_words == []:
        break
      else:
        linput_ch = int(linput_words[0].split('_')[1])
        linput_pos = int(linput_words[1])
        linput_gt = linput_words[2:]

    # introduce missing data if there is no overlap
    if pos != linput_pos:
      linput_gtNs = ['N' for l in linput_gt]
      wordsP = '\t'.join(str(e) for e in words+linput_gtNs)
    else:
      wordsP = '\t'.join(str(e) for e in words+linput_gt)

    output.write('%s\n' % wordsP)

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

linput.close()
output.close()
print('Done!')
