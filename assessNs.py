#!/usr/bin/python2

"""
This script calculates missing data (Ns) per position

# input:

Chr		pos	ind1	ind2	 ind3	ind4	ind5	ind6	ind7	ind8
chr1	1	A	G	G	A	G	G	R	A
chr1	2	Y	C	T	C	C	C	T	C

# output:

Two histograms showing the distribution on Ns according to samples and sites.

# command

$ python assessNs.py -i datafile -o outputfile -s "ind1,ind2,ind3,ind4"

contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""
############################# modules #############################

import argparse, sys, re, collections
import matplotlib
matplotlib.use('Agg') # to avoid RuntimeError('Invalid DISPLAY variable'). Must be before importing matplotlib.pyplot!
import matplotlib.pyplot as plt
import numpy as np

############################# options #############################

class MyParser(argparse.ArgumentParser): 
   def error(self, message):
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

parser = MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples for with to calculate Ns', type=str, required=True)
args = parser.parse_args()

counter = 0

############################# program #############################

print('Opening the file...')
with open(args.input) as datafile:
  header_line = datafile.readline()
  header_words = header_line.split()

  # index samples
  samplesCol = [str(j) for j in re.findall(r'[^,\s]+', args.samples)]
  numberSites = len(samplesCol)
  sampCol = []
  sampColnames = []
  sampNs = []
  for i in samplesCol:
    indnumber = header_words.index(i)
    sampCol.append(indnumber)
    sampColnames.append(header_words[indnumber])
    sampNs.append(0)

  print('Counting Ns ...')

  Ns = []

  for line in datafile:
    words = line.split()
    chr_pos = words[0:2]

    # select samples
    sample_charaters = []
    for el in sampCol:
      sample_charaters.append(words[el])

    # count Ns per position
    count=collections.Counter(sample_charaters)
    contNsOnly = count['N']
    Ns.append(contNsOnly)

    # count Ns per sample
    for i in range(len(sampCol)):
      if sample_charaters[i] == 'N':
        sampNs[i] += 1

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()

# Plot the results 
bins = np.arange(len(sampCol)) - 0.5

# Plot the histogram of Ns per position
plt.hist(Ns, color="grey", bins=bins)
plt.xlim(-1.5, len(sampCol)+1.5)
plt.xticks(range(0,len(sampCol)+1, 2))
plt.xlabel("Number of Ns")
plt.ylabel("Number of sites")
plt.title("Ns per site", size = 18)
plt.tight_layout()
plt.savefig(args.input+"_Ns_per_site.png", dpi=90)
plt.close()

# Plot the barplot of Ns per sample
plt.bar(bins, sampNs, width = 1, color = "grey", align = 'center')
plt.xticks(bins, sampColnames, rotation=90)
plt.xlim([-2,len(sampColnames)])
plt.tick_params(axis='both', labelsize=8)
plt.ylabel("Number of Ns")
plt.title("Ns per sample", size = 18)
plt.tight_layout()
plt.savefig(args.input+"_Ns_per_sample.png", dpi=90)
plt.close()

print('Done!')
