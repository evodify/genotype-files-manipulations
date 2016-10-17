#!/usr/bin/python2

""" *********************************** calls to ped *****************************************
This script converts calls file to ped file

The input is expected to be phased. So each sample except REF is presented by two sequences distinguished by _1 and _2.

The input file should look like this:
Chr	pos	REF  indA_1	indA_2	indB_1	indB_2	indC_1	 indC_2
HE669513	1	A	G	G	A	T	G	G
HE669513	2	G	C	G	C	C	C	T



The pedInfo file should look like this:

ref REF 0 0 0 0
gr1 indA 0 0 0 0
gr2 indB 0 0 0 0
gr2 indC 0 0 0 0

The output will look like this:

ref REF 0 0 0 0	A A	G G
gr1 indA 0 0 0 0	G G	C G
gr2 indB 0 0 0 0	A T	C C
gr2 indC 0 0 0 0	G G	C T

command:
$ python2 calls_to_ped.py -i <input> -o <output> -p <ped_info>

contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

import sys, argparse

########################### input options ##########################

class MyParser(argparse.ArgumentParser): 
   def error(self, message):
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

parser = MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-p', '--pedInfo', help = 'first 6 column of a ped file', type=str, required=True)
args = parser.parse_args()

############################# functions #############################

def print_GT(z, g, gg):
  gr = ' '.join(str(e) for e in groups[z])
  outFile.write("%s\t" % gr)
  for i in range(len(g)):
    if i == len(g)-1:
      outFile.write("%s %s\n" % (g[i], gg[i]))
    else:
      outFile.write("%s %s\t" % (g[i], gg[i]))
      
############################ main script ############################

ped = open(args.pedInfo, "r")

# index samples from pedInfo
pedgr = []
pedid = []
for line in ped:
  pedgr.append(line.split()[1])
  pedid.append(line.split())
groups = {}
for i in set(pedgr):
  for j in pedid:
    if i in j:
      groups[i] = j
      
# create output file
outFile = open(args.output, "w")

# convert file
with open(args.input) as datafile:
  lis = [line.split()[2:] for line in datafile] # read all lines into memory
  
  # process REF, which is haploid
  z = zip(*lis)[0][0]
  ref = zip(*lis)[0][1:]
  print_GT(z, ref, ref)
  
  # process other genotypes which are phased (2 sequence per sample)
  z1 = zip(*lis)[1][0].split("_")[0]
  g1 = zip(*lis)[1][1:]
  g2 = zip(*lis)[2][1:]

  for x in zip(*lis)[1:]:
    z2 = x[0].split("_")[0]
    if z1 != z2:
      print_GT(z1, g1, g2)
    z1 = z2[:]
    g1 = g2[:]
    g2 = x[1:]
    
  # print last string that is stored in the memory
  print_GT(z1, g1, g2)
