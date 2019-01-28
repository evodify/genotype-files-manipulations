#!/usr/bin/env python
'''
This script extracts info coded with intervals for any data with CHROM POS coordinates.

scores.file:

CHROM	startPOS	endPOS	12.4.SNPs.haplotype
scaffold_1	679	733	0.0
scaffold_1	1733	7988	0.0
scaffold_1	8339	10393	0.0
scaffold_1	10827	11638	0.0
scaffold_1	12147	14440	0.0
scaffold_1	14986	41140	0.0
scaffold_1	41583	55939	0.0
scaffold_1	56678	61321	0.602739726027
scaffold_1	61594	62590	0.0
scaffold_1	63013	63032	0.0
scaffold_1	63186	63915	0.0
scaffold_1	64331	69229	0.0
scaffold_1	69674	71548	0.0
scaffold_1	72356	73560	0.0
scaffold_1	73952	74255	1.0
scaffold_1	76025	76100	NA
scaffold_1	76499	76500	1.0
scaffold_1	77993	78223	0.0
scaffold_1	78581	78956	0.0
scaffold_1	79348	79750	NA
scaffold_1	80581	80986	NA
scaffold_1	81577	81767	0.0
scaffold_1	82218	85489	0.0
scaffold_1	85997	86249	NA
scaffold_1	87071	88130	0.0
scaffold_1	88973	89426	0.0
scaffold_1	89762	90013	0.0
scaffold_1	90564	90680	0.0
scaffold_1	91067	91347	0.0
scaffold_1	91724	100495	0.0
scaffold_1	101065	111027	0.0
scaffold_2	6186	6904	NA
scaffold_2	8480	11480	0.0
scaffold_2	11941	15733	0.0

genotype.file:

CHROM	POS	Value
scaffold_1	704	0.000000
scaffold_1	1704	0.000000
scaffold_1	2704	0.462166
scaffold_1	3704	0.616884
scaffold_1	4704	0.375514
scaffold_1	5704	0.190936
scaffold_1	6704	0.000000
scaffold_1	7704	0.000218
scaffold_1	8704	0.000246
scaffold_1	97041	0.076171
scaffold_2	1762	0.021667
scaffold_2	2762	0.405055
scaffold_2	3762	0.754771
scaffold_2	4762	0.429861
scaffold_2	5762	0.377096
scaffold_2	6762	1.022632
scaffold_2	7762	0.857520
scaffold_2	8762	0.021991
scaffold_2	9762	0.336997

output.file:

CHROM	POS	Value	score
scaffold_1	704	0.000000	0.0
scaffold_1	2704	0.462166	0.0
scaffold_1	3704	0.616884	0.0
scaffold_1	4704	0.375514	0.0
scaffold_1	5704	0.190936	0.0
scaffold_1	6704	0.000000	0.0
scaffold_1	7704	0.000218	0.0
scaffold_1	8704	0.000246	0.0
scaffold_1	97041	0.076171	0.0
scaffold_2	6762	1.022632	NA
scaffold_2	8762	0.021991	0.0
scaffold_2	9762	0.336997	0.0

command:

python select_interval_info_for_genotypes.py -i genotype.file -s scores.file -o output.file


contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

############################ options ##############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the file with genomic coordinates and values', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s', '--scores_file', help = 'phasing scores (output of assign_HapCUT_blocks.py)', type=str, required=True)
args = parser.parse_args()

############################# program #############################

scoresFile = open(args.scores_file, "r")
scoresWords = scoresFile.readline().split()
scoresScore = scoresWords[3:]
scoresWords = scoresFile.readline().split()
scoresScaf = scoresWords[0]
scoresStart = int(scoresWords[1])
scoresEnd = int(scoresWords[2])

output = open(args.output, 'w')

counter = 0

with open(args.input) as input:
  header_line = input.readline()
  header_words = header_line.split()
  headerP = '\t'.join(str(e) for e in header_words)
  scoresScoreP = '\t'.join(str(e) for e in scoresScore)

  output.write("%s\t%s\n" % (headerP, scoresScoreP))

  for line in input:
    inputWords = line.split()
    inputScaf = inputWords[0]
    inputPos = int(inputWords[1])
    inputVal = str(inputWords[2])

    while ((inputScaf != scoresScaf) or (inputScaf == scoresScaf and inputPos > scoresEnd)):
      scoresWords = scoresFile.readline().split()
      if not scoresWords:
        break
      else:
        scoresScaf = scoresWords[0]
        scoresStart = int(scoresWords[1])
        scoresEnd = int(scoresWords[2])
        scoresScore = scoresWords[3:]


    if (inputScaf == scoresScaf and inputPos >= scoresStart and inputPos <= scoresEnd):
      inputWordsP = '\t'.join(str(e) for e in inputWords)
      if len(scoresScore) > 1:
        scoresScoreP = '\t'.join(str(e) for e in scoresScore)
      else:
        scoresScoreP = scoresScore
      output.write('%s\t%s\n' % (inputWordsP, scoresScoreP))

    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

input.close()
scoresFile.close()
output.close()
print('Done!')
