#! /usr/bin/env python
"""
This script merges results of a sliding window analysis scores with every SNPs.


# SNPs.txt:

chr	POS	
6	16
6	43
6	66
6	68
6	94


# scores.txt:

chr	POS	score
6	10	0.0010
6	20	0.0021
6	30	0.0031
6	40	0.0041
6	50	0.0051
6	60	0.0061
6	70	0.0072
6	80	0.0082
6	90	0.0092

# output.txt:

chr	POS	score
chr6	16	0.0021
chr6	43	0.0041
chr6	66	0.0072
chr6	68	0.0072
chr6	94	0.0092

# command:

$ python merge_windowScores_SNPs.py \
    -i SNPs.txt \
    -s scores.txt \
    -o output.tab

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""
############################# modules #############################

import calls  # my custom module

############################# options #############################

parser = calls.CommandLineParser()
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

counter = 0

roundNumber = -2

scoresFile = open(args.scores, 'r')
scoresFile_header = scoresFile.readline().split()
scoresFile_headerP = '\t'.join(str(e) for e in scoresFile_header[2:])
scoresFile_words = scoresFile.readline().split()
scoresFile_CHR = int(scoresFile_words[0].split('chr')[-1])
scoresFile_POS = int(round(float(scoresFile_words[1]), roundNumber))
scoresFile_content = scoresFile_words[2:]

print('Opening the file...')
with open(args.input) as datafile:
    header_line = datafile.readline().split()
    header_lineP = '\t'.join(str(e) for e in header_line[0:2])
    output = open(args.output, 'w')
    output.write("%s\t%s\n" % (header_lineP, scoresFile_headerP))

    for line in datafile:
        words = line.split()
        ch = int(words[0].split('chr')[-1])
        pos = int(words[1])
        posR = int(round(float(words[1]), roundNumber))

        pCHR = scoresFile_CHR
        pPOS = scoresFile_POS

        # read the score file if necessary
        while ((ch == scoresFile_CHR and posR > scoresFile_POS) or \
                (ch > scoresFile_CHR)) and (scoresFile_words != []):
                scoresFile_words = scoresFile.readline().split()
                if scoresFile_words != []:    
                    scoresFile_CHR = int(scoresFile_words[0].split('chr')[-1])
                    scoresFile_POS = int(round(float(scoresFile_words[1]), roundNumber))
                    scoresFile_content = scoresFile_words[2:]

        # find overlap
        if ch == scoresFile_CHR and posR == scoresFile_POS:
            contentP =  '\t'.join(str(e) for e in scoresFile_content)
            output.write('chr%s\t%s\t%s\n' % (ch, pos, contentP))
        else:
            print "WARNING: No match found for chr%s %s"  %  (ch, pos)
            print ch, pos, posR, ":",  scoresFile_CHR, scoresFile_POS, "=", pCHR, pPOS
  
        counter = calls.lineCounter(counter)

datafile.close()
scoresFile.close()
output.close()

print('Done!')
