#!/usr/bin/env python
'''
This script extracts gene names from a bed file by provided coordinates.

bed.file:


interval.file:


output.file:


command:

python select_genes_by_intervals.py -i interval.file -b bed.file -o output.file


contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--interval_file', help = 'name of the file with genome intervals', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-b', '--bed_file', help = 'file containing list of genes with scaffolds and position information', type=str, required=True)
args = parser.parse_args()

############################# program #############################

interFile = open(args.interval_file, "r")
interHeader = interFile.readline()
interHeaderP = '\t'.join(str(e) for e in interHeader.split())
IntervalWords = interFile.readline().split()
IntervalScaf = int(IntervalWords[0].split('_')[1])
IntervalStart = int(IntervalWords[1])
IntervalEnd = int(IntervalWords[2])

output = open(args.output, 'w')
output.write("%s\tgenes\n" % interHeaderP)

counter = 0

geneL = []

with open(args.bed_file) as bedfile:
  for line in bedfile:
    geneWords = line.split()
    geneScaf = int(geneWords[0].split('_')[1])
    geneStart = int(geneWords[1])
    geneEnd = int(geneWords[2])
    geneName = str(geneWords[3])


    if (geneScaf == IntervalScaf) and ((IntervalStart <= geneStart and geneEnd <= IntervalEnd) or (IntervalStart <= geneStart and geneStart <= IntervalEnd) or (IntervalStart <= geneEnd and geneEnd <= IntervalEnd)):
      geneL.append(geneName)
    elif (IntervalScaf < geneScaf) or (IntervalScaf == geneScaf and IntervalStart < geneStart):
      IntervalWordsP = '\t'.join(str(e) for e in IntervalWords)
      if geneL==[]:
        geneLP = "NA"
      else:
        geneLP = ','.join(str(e) for e in geneL)
      output.write('%s\t%s\n' % (IntervalWordsP, geneLP))
      geneL = []
      IntervalWords = interFile.readline().split()
      if IntervalWords == []:
        break
      else:
        IntervalScaf = int(IntervalWords[0].split('_')[1])
        IntervalStart = int(IntervalWords[1])
        IntervalEnd = int(IntervalWords[2])
      if (geneScaf == IntervalScaf) and ((IntervalStart <= geneStart and geneEnd <= IntervalEnd) or (IntervalStart <= geneStart and geneStart <= IntervalEnd) or (IntervalStart <= geneEnd and geneEnd <= IntervalEnd)):
        geneL.append(geneName)
    #else: # for debugging
      #print "interval:", IntervalScaf, IntervalStart, IntervalEnd, "BED:", geneScaf, geneStart, geneEnd

    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

if IntervalWords != []:
  IntervalWordsP = '\t'.join(str(e) for e in IntervalWords)
  if geneL==[]:
    geneLP = "NA"
  else:
    geneLP = ','.join(str(e) for e in geneL)
  output.write('%s\t%s\n' % (IntervalWordsP, geneLP))

interFile.close()
bedfile.close()
output.close()
print('Done!')
