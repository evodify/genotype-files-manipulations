#! /usr/bin/env python2

'''
Annotates genes from a sliding windows analysis with the stats per gene.
It handles overlapping intervals and genes spanning a few windows by outputting mean, max and min stats.


input.txt:

chr	start	end	stats1	stats2	genes
chr1	1	10	0.527	5.3	ENSCAFG00000000170,ENSCAFG00000000171,ENSCAFG00000000172,ENSCAFG00000029562,ENSCAFG00000029833
chr1	2	20	0.364	3.6	ENSCAFG00000000171,ENSCAFG00000000172,ENSCAFG00000029562,ENSCAFG00000029833,ENSCAFG00000000173,ENSCAFG00000031536,ENSCAFG00000023012
chr1	3	30	0.501	5	ENSCAFG00000000239
chr1	4	40	0.405	4.1	ENSCAFG00000000239
chr2	10	20	0.402	4	ENSCAFG00000031480,ENSCAFG00000000894,ENSCAFG00000000897,ENSCAFG00000000901,ENSCAFG00000000905,ENSCAFG00000000908,ENSCAFG00000000911
chr2	20	30	0.443	4.4	ENSCAFG00000001026
chr2	30	40	0.453	4.5	ENSCAFG00000001026
chr2	40	50	0.43	4.3	ENSCAFG00000001026
chr2	50	60	0.38	3.8	ENSCAFG00000030357

output.txt:

gene	mean_stats1	max_stats1	min_stats1	mean_stats2	max_stats2	min_stats2
ENSCAFG00000029562	0.4455	0.527	0.364	4.45	5.3	3.6
ENSCAFG00000000911	0.402	0.402	0.402	4.0	4.0	4.0
ENSCAFG00000000897	0.402	0.402	0.402	4.0	4.0	4.0
ENSCAFG00000000894	0.402	0.402	0.402	4.0	4.0	4.0
ENSCAFG00000029833	0.4455	0.527	0.364	4.45	5.3	3.6
ENSCAFG00000023012	0.364	0.364	0.364	3.6	3.6	3.6
ENSCAFG00000030357	0.38	0.38	0.38	3.8	3.8	3.8
ENSCAFG00000031536	0.364	0.364	0.364	3.6	3.6	3.6
ENSCAFG00000000908	0.402	0.402	0.402	4.0	4.0	4.0
ENSCAFG00000001026	0.442	0.453	0.43	4.3999999999999995	4.5	4.3
ENSCAFG00000000905	0.402	0.402	0.402	4.0	4.0	4.0
ENSCAFG00000031480	0.402	0.402	0.402	4.0	4.0	4.0
ENSCAFG00000000901	0.402	0.402	0.402	4.0	4.0	4.0
ENSCAFG00000000239	0.453	0.501	0.405	4.55	5.0	4.1
ENSCAFG00000000171	0.4455	0.527	0.364	4.45	5.3	3.6
ENSCAFG00000000170	0.527	0.527	0.527	5.3	5.3	5.3
ENSCAFG00000000173	0.364	0.364	0.364	3.6	3.6	3.6
ENSCAFG00000000172	0.4455	0.527	0.364	4.45	5.3	3.6

#command:

$ python2 annotate_genes_withSlidingWindowsStats.py \
    -i input.txt \
    -o output.txt

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu
'''

############################# modules #############################

import calls # my custom module
import numpy as np

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
args = parser.parse_args()

############################# program #############################

outfile = open(args.output, 'w')

with open(args.input) as datafile:
  header_words = datafile.readline().split()
  statName = header_words[3:-1]

  output = open(args.output, 'w')
  output.write('gene')
  for i in statName:
    output.write('\tmean_%s\tmax_%s\tmin_%s' % (i,i,i))
  output.write('\n')

  genesDics = {}

  for line in datafile:
    words = line.split()
    stats = []
    for i in words[3:-1]:
      if i=='NA':
        stats.append(np.nan)
      else:
        stats.append(float(i))
    genes = words[-1].split(',')
    for g in genes:
      if g in genesDics.keys():
        genesDics[g].append(stats)
      else:
        genesDics[g] = [stats]
  for k in genesDics.keys():
    mt = np.stack(genesDics[k])
    mt_mean = np.nanmean(mt, axis=0)
    mt_min = np.nanmin(mt, axis=0)
    mt_max = np.nanmax(mt, axis=0)
    
    output.write(k)
    for i in range(0,len(statName)):
      output.write('\t%s\t%s\t%s' % (mt_mean[i],mt_max[i],mt_min[i]))
    output.write('\n')

datafile.close()
outfile.close()

print('Done!')
