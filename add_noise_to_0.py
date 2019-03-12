#!/usr/bin/python2
"""
This script replaces 0 with very small noise values.

# input:

chr	POS	score
chr6	16	0.0000
chr6	43	0.0041
chr6	66	0.0000
chr6	68	0.0072
chr6	94	0.0092


# output:

chr	POS	score
chr6	16	3.05121406663e-06
chr6	43	0.0041
chr6	66	-2.97204390375e-06
chr6	68	0.0072
chr6	94	0.0092

# command

$ python add_noise_to_0.py -i input -o outputfile -s "score"

contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

############################# modules #############################

import numpy as np
import calls  # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument(
    '-i',
    '--input',
    help='name of the input file',
    type=str,
    required=True)
parser.add_argument(
    '-o',
    '--output',
    help='name of the output file',
    type=str,
    required=True)
parser.add_argument(
    '-s',
    '--stats',
    help='column names of the stats to process',
    type=str,
    required=True)
args = parser.parse_args()

sampleNames = calls.checkSampleNames(args.stats, args.input)

############################# program #############################

output = open(args.output, 'w')

counter = 0

print('Opening the file...')

with open(args.input) as datafile:
    header_line = datafile.readline()
    header_words = header_line.split()
    
    sampCol = calls.indexSamples(sampleNames, header_words)
    sampColnames = calls.selectSamples(sampCol, header_words)
    sampColnamesP =  '\t'.join(str(e) for e in sampColnames)
    
    output.write("%s\t%s\t%s\n" %
                (header_words[0], header_words[1], sampColnamesP))

    print('Adding noise ...')

    for line in datafile:
        words = line.split()
        chr_pos = words[0:2]
        sample_scores = calls.selectSamples(sampCol, words)

        sample_scoresNoise = []

        for s in sample_scores:
            if float(s) == 0.0:
                noise = np.random.normal(0.000001,0.00001,1)[0]
                sample_scoresNoise.append(noise)
            else:
                sample_scoresNoise.append(s)

        sample_scoresNoiseP =  '\t'.join(str(e) for e in sample_scoresNoise)
        output.write("%s\t%s\t%s\n" %
                (words[0], words[1], sample_scoresNoiseP))
        
        counter = calls.lineCounter(counter)

datafile.close()
output.close()
print('Done!')
