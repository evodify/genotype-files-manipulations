#! /usr/bin/env python2

'''
This script annotates genes from a sliding windows analysis with the stats per gene.
It handles overlapping intervals and genes spanning a few windows by outputting mean, max and min stats.

#Example input:

stats                   genes
0.527270908085201       ENSCAFG00000000170,ENSCAFG00000000171,ENSCAFG00000000172,ENSCAFG00000029562,ENSCAFG00000029833
0.364324009629718       ENSCAFG00000000171,ENSCAFG00000000172,ENSCAFG00000029562,ENSCAFG00000029833,ENSCAFG00000000173,ENSCAFG00000031536,ENSCAFG00000023012
0.500663017366094       ENSCAFG00000000239
0.405364585439997       ENSCAFG00000000239
0.402080897939156       ENSCAFG00000031480,ENSCAFG00000000894,ENSCAFG00000000897,ENSCAFG00000000901,ENSCAFG00000000905,ENSCAFG00000000908,ENSCAFG00000000911
0.443099011567851       ENSCAFG00000001026
0.452831729053616       ENSCAFG00000001026
0.430122661190792       ENSCAFG00000001026
0.379507893057865       ENSCAFG00000030357


#Example output:

genes	mean_stats	max_stats	min_stats
ENSCAFG00000029562	0.445797458857	0.527270908085	0.36432400963
ENSCAFG00000000911	0.402080897939	0.402080897939	0.402080897939
ENSCAFG00000000897	0.402080897939	0.402080897939	0.402080897939
ENSCAFG00000000894	0.402080897939	0.402080897939	0.402080897939
ENSCAFG00000029833	0.445797458857	0.527270908085	0.36432400963
ENSCAFG00000023012	0.36432400963	0.36432400963	0.36432400963
ENSCAFG00000030357	0.379507893058	0.379507893058	0.379507893058
ENSCAFG00000031536	0.36432400963	0.36432400963	0.36432400963
ENSCAFG00000000908	0.402080897939	0.402080897939	0.402080897939
ENSCAFG00000001026	0.442017800604	0.452831729054	0.430122661191
ENSCAFG00000000905	0.402080897939	0.402080897939	0.402080897939
ENSCAFG00000031480	0.402080897939	0.402080897939	0.402080897939
ENSCAFG00000000901	0.402080897939	0.402080897939	0.402080897939
ENSCAFG00000000239	0.453013801403	0.500663017366	0.40536458544
ENSCAFG00000000171	0.445797458857	0.527270908085	0.36432400963
ENSCAFG00000000170	0.527270908085	0.527270908085	0.527270908085
ENSCAFG00000000173	0.36432400963	0.36432400963	0.36432400963
ENSCAFG00000000172	0.445797458857	0.527270908085	0.36432400963


#command:

$ python annotate_genes_withSlidingWindowsStats.py -i input.table -o output.tab

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
args = parser.parse_args()

############################# program #############################

outfile = open(args.output, 'w')

with open(args.input) as datafile:
  header_words = datafile.readline().split()
  output = open(args.output, 'w')
  output.write('genes\tmean_%s\tmax_%s\tmin_%s\n' % (header_words[0], header_words[0], header_words[0]))
  genesDics = {}

  for line in datafile:
    words = line.split()
    stats = words[0]
    genes = words[1].split(',')
    for g in genes:
      if g in genesDics.keys():
        genesDics[g].append(stats)
      else:
        genesDics[g] = [stats]
  for k in genesDics.keys():
    v = []
    for i in genesDics[k]:
      if i != "NA":
        v.append(float(i))
    try:
      output.write('%s\t%s\t%s\t%s\n' % (k, sum(v)/len(v), max(v), min(v)))
    except:
      print("WARNING: Stats for %s is not numeric. The output value will be NA" % k)
      output.write('%s\tNA\tNA\n' % k)

datafile.close()
outfile.close()

print('Done!')
