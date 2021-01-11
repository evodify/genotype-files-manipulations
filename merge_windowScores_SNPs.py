#! /usr/bin/env python
"""
This script merges results of a sliding window analysis scores with every SNPs.


# SNPs.txt:

chr	POS	
chr6	16
chr6	43
chr6	66
chr6	68
chr6	94


# scores.txt:

chr	POS	score
chr6	10	0.0010
chr6	20	0.0021
chr6	30	0.0031
chr6	40	0.0041
chr6	50	0.0051
chr6	60	0.0061
chr6	70	0.0072
chr6	80	0.0082
chr6	90	0.0092

# output.txt:

chr	POS	score
chr6	16	0.0021
chr6	43	0.0041
chr6	66	0.0072
chr6	68	0.0072
chr6	94	0.0092

# command:

$ python2 merge_windowScores_SNPs.py \
    -i SNPs.txt \
    -s scores.txt \
    -o output.tab

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""
############################# modules #############################

import argparse, sys  # for input options

############################# options #############################

class CommandLineParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = CommandLineParser()
parser.add_argument(
    '-i',
    '--input',
    help='name of the input SNPs file',
    type=str,
    required=True)
parser.add_argument(
    '-s',
    '--scores',
    help='name of the input scores files',
    type=str,
    required=True)
parser.add_argument(
    '-o',
    '--output',
    help='name of the output file',
    type=str, required=True)
args = parser.parse_args()

############################# program #############################

roundNumber = 'NA'
n = 1

output = open(args.output, 'w')

with open(args.scores) as scoresFile: 
  scoresHeader = scoresFile.readline()
  output.write(scoresHeader)
  scoresDic = {}
  for line in scoresFile:
    scoresWords = line.rstrip().split()
    if roundNumber == 'NA' and n == 1:
        firstScoreChr = scoresWords[0]
        firstScorePos = int(scoresWords[1])
        firstScores = scoresWords[2:]
        n+=1
    elif roundNumber == 'NA' and n == 2:
        roundNumber = -(len(str(int(scoresWords[1])-firstScorePos))-1)
        scoresPOS = int(round(float(scoresWords[1]), roundNumber))
        scoresCoord = str(scoresWords[0]) + ":" + str(scoresPOS)
        scores = scoresWords[2:]
        scoresDic[scoresCoord] = scores

        firstScorePosRound = int(round(float(firstScorePos), roundNumber))
        firstScoreCoord = str(firstScoreChr) + ":" + str(firstScorePosRound)
        scoresDic[firstScoreCoord] = firstScores
    else:
        scoresPOS = int(round(float(scoresWords[1]), roundNumber))
        scoresCoord = str(scoresWords[0]) + ":" + str(scoresPOS)
        scores = scoresWords[2:]
        scoresDic[scoresCoord] = scores

with open(args.input) as chrPos:
    chrPosHeader = chrPos.readline()
    for line in chrPos:
      words = line.rstrip().split('\t')
      scoresPOS = int(round(float(words[1]), roundNumber))
      chrPosCoord = str(words[0]) + ":" + str(scoresPOS)
      try:
          annotScores = scoresDic[chrPosCoord][:]
      except:
          annotScores = ['NA']*len(firstScores)
      annotScoresP = '\t'.join(str(e) for e in annotScores)
      output.write("%s\t%s\n" % (line.rstrip(), annotScoresP))

print('Done!')
