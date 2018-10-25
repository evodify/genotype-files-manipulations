#!/usr/bin/env python

'''
This script merges all chromosomes into continuous genomic coordinates.

Example input:

CHROM   POS SomeValues
chr_1   1   2456
chr_1   2   36
chr_1   3   346
chr_1   4   36
chr_2   1   36
chr_2   2   2461
chr_2   3   6
chr_2   4   70
chr_2   5   2464
chr_3   1   46
chr_3   2   2466
chr_3   3   36
chr_3   4   6
chr_3   5   6

Example fasta.fai (sorted by coordinates):

chr_1      19624517
chr_2      14106692
chr_3      15060676

Example output:

Genome_Pos	CHROM	POS	SomeValues
1.0	chr_1	1	2456
2.0	chr_1	2	36
3.0	chr_1	3	346
4.0	chr_1	4	36
19624518.0	chr_2	1	36
19624519.0	chr_2	2	2461
19624520.0	chr_2	3	6
19624521.0	chr_2	4	70
19624522.0	chr_2	5	2464
33731210.0	chr_3	1	46
33731211.0	chr_3	2	2466
33731212.0	chr_3	3	36
33731213.0	chr_3	4	6
33731214.0	chr_3	5	6

command:

$ python mergeChrPos_in_callsTab.py -i input.tab -f fasta.fai -o output.file

contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''
############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-f', '--fai', help = 'name of the fasta.fai file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
args = parser.parse_args()

############################# program #############################

counter = 0

print('Opening files...')

# create chromosome whole genomes position coordinates
fai = open(args.fai, 'r')
faiDict = {}
gPosition = 0
for line in fai:
  fai_words = line.split()
  faiDict[fai_words[0]] =  gPosition
  gPosition += int(fai_words[1])
fai.close()

output = open(args.output, 'w')

print('Opening the file...')
with open(args.input) as datafile:
  header = datafile.readline().split()
  headerP = '\t'.join(str(e) for e in header)
  output.write("Genome_Pos\t%s\n" % headerP)
  for line in datafile:
    words = line.split()
    chr = words[0]
    pos = float(words[1])
    wordsP = '\t'.join(str(e) for e in words)

    if chr not in faiDict.keys():  # end of fai file
      raise Exception("%s is not found in %s.\nExecution stopped!" % (chr, args.fai))
    else:
      fai_start = faiDict[chr]
    posP = pos+fai_start

    output.write("%s\t%s\n" % (posP, wordsP))

    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
output.close()
print('Done!')
 
