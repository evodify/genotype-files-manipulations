#! /usr/bin/env python

"""
This script removes all sites that consists of more than given amount of missing data (Ns).

# input:

Chr     pos ind1    ind2     ind3   ind4    ind5    ind6    ind7    ind8
chr1    1   A   G   G   A   G   G   R   A
chr1    2   Y   C   N   C   C   C   T   C
chr1    3   Y   C   T   C   C   C   T   C
chr1    4   N   N   N   N   C   C   T   C
chr1    5   G   G   T   G   G   G   T   G

# output:
Chr     pos ind1    ind2     ind3   ind4    ind5    ind6    ind7    ind8
chr1    1   A   G   G   A   G   G   R   A
chr1    2   Y   C   N   C   C   C   T   C
chr1    3   Y   C   T   C   C   C   T   C
chr1    5   G   G   T   G   G   G   T   G

# command:

$ python filterByNs.py -i datafile -o outputfile -m 4 -s "ind1,ind2,ind3,ind4,ind5,ind6,ind7,ind8"

In this example, sites with 4 or more Ns will be removed.

contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""
############################# modules #############################

import sys, argparse, collections
import calls # my custom module

############################# options #############################

class MyParser(argparse.ArgumentParser): 
   def error(self, message):
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

parser = MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-m', '--missing', help = 'missing data threshold to remove sites', type=int, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process', type=str, required=True)

args = parser.parse_args()

############################# program #############################

counter = 0

print('Opening the file...')

with open(args.input) as datafile:
  header_line = datafile.readline()
  header_words = header_line.split()
  
  # index samples
  sampCol = calls.indexSamples(args.samples, header_words)
  
  print('Creating the output file...')
  fileoutput = open(args.output, 'w')
  fileoutput.write(header_line)
  
  for line in datafile:
    words = line.split()
    chr_pos = words[0:2]
    
    # select samples
    genotypes = calls.selectSamples(sampCol, words)

    count=collections.Counter(genotypes)
    valueN = count['N']
    
    if valueN < args.missing:
      fileoutput.write(line)

    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
fileoutput.close()

print('Done!')
