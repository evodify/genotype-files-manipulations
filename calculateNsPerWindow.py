#! /usr/bin/env python

'''
This script calculates number of positions with missing data (Ns) using sliding window approach.

#Example input:

CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   1   A   W   N   N   A   N   N   N   N
chr_1   2   C   Y   Y   N   C   C   N   C   N
chr_1   3   C   N   C   N   C   C   C   C   C
chr_1   4   T   T   T   N   T   T   T   T   T
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

#Example output:

CHROM   POS Ns
chr_1   2.5 3.0
chr_1   6.0 2.0
chr_2   3.0 2.0
chr_3   3.0 4.2
chr_4   1.5 2.0

#command:

$ python calculateNsPerWindow.py -i input.tab -o output.tab -w 5 -s "sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8"

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import sys, argparse # for input options
import calls # my custom module

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process', type=str, required=True)
parser.add_argument('-w', '--window', help = 'sliding window size', type=int, required=True)
args = parser.parse_args()

############################ functions ###########################

def processWindow(Chr, FirstPos, LastPos, Ns, outputFile):
  posP = float(FirstPos)+((float(LastPos)-float(FirstPos))/2.0)
  NsInWindow = sum(Ns)/len(Ns)
  outputFile.write("%s\t%s\t%s\n" % (Chr, posP, NsInWindow))

############################# program #############################

print('Opening the file...')

windSize = args.window
counter = 0

with open(args.input) as datafile:
  header_line = datafile.readline()
  header_words = header_line.split()

  # index samples
  sampCol = calls.indexSamples(args.samples, header_words)

  # make output header
  print('Creating the output file...')
  NsFr_output = open(args.output, 'w')
  NsFr_output.write("CHROM\tPOS\tNs\n")

############################## perform counting ####################

  print('Counting Ns ...')

  sites = 0
  Nwindow = []
  Chr2 = ''
  posS = ''
  posE = ''
  for line in datafile:
    words = line.split()
    Chr1 = words[0]
    pos = words[1]

    # to store the values of a previous line
    if not Chr2:
      Chr2 = Chr1
    if not posS:
      posS = pos
    if not posE:
      posE = pos

    # select samples
    sample_charaters = calls.selectSamples(sampCol, words)

    # if window size is reached output the results
    sites += 1
    if Chr1 != Chr2:
      processWindow(Chr2, posS, posE, Nwindow, NsFr_output)
      sites = 0
      Nwindow = []
      posS = pos
    elif sites == windSize:
      processWindow(Chr1, posS, posE, Nwindow, NsFr_output)
      sites = 0
      Nwindow = []
      posS = pos
    Chr2 = Chr1
    posE = pos

    # count Ns
    contNsOnly = calls.countPerPosition(sample_charaters, 'N')
    Nwindow.append(float(contNsOnly))

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

# process the last window
processWindow(Chr1, posS, posE, Nwindow, NsFr_output)

datafile.close()
NsFr_output.close()
print('Done!')

