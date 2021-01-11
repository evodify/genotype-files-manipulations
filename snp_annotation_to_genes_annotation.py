#! /usr/bin/env python2
'''
Transforms stats annotated with genes to genes annotated with stats.

input.txt:

chr	pos	annotation1	annotation2	genes
chr1	1	liver	liverTF	ENSCAFG00000000170,ENSCAFG00000000171,ENSCAFG00000000172,ENSCAFG00000029562,ENSCAFG00000029833
chr1	2	liver	pancreaseTF	ENSCAFG00000000171,ENSCAFG00000000172,ENSCAFG00000029562,ENSCAFG00000029833,ENSCAFG00000000173,ENSCAFG00000031536,ENSCAFG00000023012
chr1	3	pancrease	pancreaseTF	ENSCAFG00000000239
chr1	4	muscle	pancreaseTF	ENSCAFG00000000239
chr2	10	pancrease	pancreaseTF	ENSCAFG00000031480,ENSCAFG00000000894,ENSCAFG00000000897,ENSCAFG00000000901,ENSCAFG00000000905,ENSCAFG00000000908,ENSCAFG00000000911
chr2	11	muscle	muscleTF	ENSCAFG00000001026
chr2	12	muscle	muscleTF	ENSCAFG00000001026
chr2	13	muscle	muscleTF	ENSCAFG00000001026
chr2	14	muscle	muscleTF	ENSCAFG00000030357


output.txt:

gene	annotation1	annotation2
ENSCAFG00000029562	liver	pancreaseTF,liverTF
ENSCAFG00000000911	pancrease	pancreaseTF
ENSCAFG00000000897	pancrease	pancreaseTF
ENSCAFG00000000894	pancrease	pancreaseTF
ENSCAFG00000029833	liver	pancreaseTF,liverTF
ENSCAFG00000023012	liver	pancreaseTF
ENSCAFG00000030357	muscle	muscleTF
ENSCAFG00000031536	liver	pancreaseTF
ENSCAFG00000000908	pancrease	pancreaseTF
ENSCAFG00000001026	muscle	muscleTF
ENSCAFG00000000905	pancrease	pancreaseTF
ENSCAFG00000031480	pancrease	pancreaseTF
ENSCAFG00000000901	pancrease	pancreaseTF
ENSCAFG00000000239	pancrease,muscle	pancreaseTF
ENSCAFG00000000171	liver	pancreaseTF,liverTF
ENSCAFG00000000170	liver	liverTF
ENSCAFG00000000173	liver	pancreaseTF
ENSCAFG00000000172	liver	pancreaseTF,liverTF


#command:

$ python snp_annotation_to_genes_annotation.py \
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
  statName = header_words[2:-1]

  output = open(args.output, 'w')
  output.write('gene')
  for i in statName:
    output.write('\t%s' % i)
  output.write('\n')

  genesDics = {}

  for line in datafile:
    words = line.split()
    genes = words[-1].split(',')
    stats = []
    for i in words[2:-1]:
        stats.append(i)
    for g in genes:
      if g in genesDics.keys():
        genesDics[g].append(stats)
      else:
        genesDics[g] = [stats]
  for k in genesDics.keys():
    mt = np.stack(genesDics[k])
    output.write(k)
    for i in range(len(statName)):
      annot_unique = set(list([row[i] for row in mt]))
      annotP = ','.join(str(e) for e in annot_unique)
      output.write('\t%s' % annotP)
    output.write('\n')

datafile.close()
outfile.close()

print('Done!')
