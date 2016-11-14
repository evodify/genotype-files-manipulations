#!/usr/bin/env python2

""" *********************************** calls to map *****************************************
This script makes a map file from a calls file.
It is assumed the genetic distance (column 3) is not known. The column 3 will be filled with 0.

The input file should look like this:

Chr	pos	REF	indA	indB	indC	 indC
chr_1	1	A	G	A	R
chr_1	2	G	C	C	T

Note! Chromosome number must be separated by _. For example, chr_1 - correct, chr1 - incorrect.
Genotypes field actually do not matter for this conversion.

The output will look like this:

1	chr_1_1	0	1
1	chr_1_2	0	2

command:
$ python2 calls_to_map.py -i <input> -o <output>

contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

import argparse, sys

class MyParser(argparse.ArgumentParser): 
   def error(self, message):
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

parser = MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
args = parser.parse_args()

output = open(args.output, 'w')

with open(args.input) as datafile:
  header = datafile.readline()
    
  for line in datafile:
    words = line.split()
    ch = int(words[0].split('_')[1])
    pos = int(words[1])
    indx = str(words[0])+'_'+str(pos)
    
    output.write("%s\t%s\t%s\t%s\n" % (ch, indx, 0, pos))
    
output.close()
