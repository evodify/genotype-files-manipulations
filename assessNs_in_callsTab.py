#!/usr/bin/python2

"""
This script calculates missing data (Ns) per position/sample and visualizes the results.

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

import matplotlib
matplotlib.use('Agg') # to avoid RuntimeError('Invalid DISPLAY variable'). Must be before importing matplotlib.pyplot!
import matplotlib.pyplot as plt
import numpy as np
import collections
import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process (optional)', type=str, required=False)
args = parser.parse_args()

# check if samples names are given and if all sample names are present in a header
sampleNames = calls.checkSampleNames(args.samples, args.input)

############################# program #############################

counter = 0

print('Opening the file...')

with open(args.input) as datafile:
  header_line = datafile.readline()
  header_words = header_line.split()

  # index samples
  sampCol = calls.indexSamples(sampleNames, header_words)

  # create lists for output
  sampColnames = calls.selectSamples(sampCol, header_words)
  sampNs = [0 for i in sampColnames]

  print('Counting Ns ...')

  Ns = []

  for line in datafile:
    words = line.split()
    chr_pos = words[0:2]

    # select samples
    sample_charaters = calls.selectSamples(sampCol, words)

    # count Ns per position
    contNsOnly = calls.countPerPosition(sample_charaters, 'N')
    Ns.append(contNsOnly)

    # count Ns per sample
    calls.countPerSample(sample_charaters, sampNs, 'N')

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()

# write the counts to a fine
outputTXTsite = open(args.output+"_Ns_per_site.csv", 'w')
outputTXTsample = open(args.output+"_Ns_per_sample.csv", 'w')

outputTXTsite.write("Number_of_Ns\tNumber_of_sites\n")
outputTXTsample.write("Sample\tNumber_of_Ns\n")

Ns.sort()
counter=collections.Counter(Ns)

for k,v in zip(counter.keys(), counter.values()):
    outputTXTsite.write("%s\t%s\n" % (k, v))

for s,n in zip(sampColnames, sampNs):
    outputTXTsample.write("%s\t%s\n" % (s, n))

# Plot the results 
print('Printing graphics ...')

binsPos = np.arange(len(sampCol)+2) - 0.5
binsSamp = np.arange(len(sampCol)) - 0.5
heightF = 6

# increase figure width if the number of samples is too large
if len(sampCol) > 50:
  widthF = float(len(sampCol))/5.0
else:
  widthF = 8

# Plot the histogram of Ns per position
plt.figure(figsize=(widthF, heightF))
plt.hist(Ns, color="grey", bins=binsPos)
plt.xlim(-1.5, len(sampCol)+1.5)
plt.xticks(range(0,len(sampCol)+1, 1))
plt.xlabel("Number of Ns")
plt.ylabel("Number of sites")
plt.title("Ns per site", size = 18)
plt.tight_layout()
plt.savefig(args.output+"_Ns_per_site.png", dpi=90)
plt.close()

# Plot the barplot of Ns per sample
plt.figure(figsize=(widthF, heightF))
plt.bar(binsSamp, sampNs, width = 1, color = "grey", align = 'center')
plt.xticks(binsSamp, sampColnames, rotation=90)
plt.xlim([-2,len(sampColnames)])
plt.tick_params(axis='both', labelsize=8)
plt.ylabel("Number of Ns")
plt.title("Ns per sample", size = 18)
plt.tight_layout()
plt.savefig(args.output+"_Ns_per_sample.png", dpi=90)
plt.close()

print('Done!')
