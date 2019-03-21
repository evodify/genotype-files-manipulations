#! /usr/bin/env python
"""
This script merges tab files by CHR and POS with replacement of the overlap by 
the input from the second input file. 

I used 10bp sliding window to get the input files, so in this script I 
round each position for the nearest 10th position. If your input differs,
modify or remove the parts of the code containing "round("

# test1.txt:

chr	POS	score	max_s
6	11	1	0.0001
6	21	2	0.0001
6	31	3	0.0001
6	41	4	0.0001
6	51	5	0.0001
6	61	6	0.0001
6	71	7	0.0001
6	81	8	0.0001
6	91	9	0.0001


# test2.txt

chr	POS	score	max_s
6	72	10	0.0002
6	82	11	0.0002
6	92	12	0.0002
6	102	13	0.0002
6	112	14	0.0002
6	122	15	0.0002
6	132	16	0.0002
6	142	17	0.0002
6	152	18	0.0002


# output:

chr	POS	score	max_s
6	10	1	0.0001
6	20	2	0.0001
6	30	3	0.0001
6	40	4	0.0001
6	50	5	0.0001
6	60	6	0.0001
6	70	10	0.0002
6	80	11	0.0002
6	90	12	0.0002


# command:

$ python merge_overlaping_windows_replaceOverlap.py \
    -r test1.txt \
    -i test2.txt \
    -o output.tab

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""
############################# modules #############################

import calls  # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument(
    '-r',
    '--reference_input',
    help='name of the reference input file',
    type=str,
    required=True)
parser.add_argument(
    '-i',
    '--input_to_merge',
    help='name of the input file to introduce into the reference input',
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

inputFile2 = open(args.input_to_merge, 'r')
inputFile2_header = inputFile2.readline().split()
inputFile2_words = inputFile2.readline().split()
inputFile2_CHR = int(inputFile2_words[0].split('chr')[-1])
inputFile2_POS = int(round(float(inputFile2_words[1]), -1))
inputFile2_content = inputFile2_words[2:]
    
output = open(args.output, 'w')

print('Opening the file...')
with open(args.reference_input) as datafile:
    header_line = datafile.readline().split()

    if header_line != inputFile2_header:
        raise IOError('Headers are different in the two files')
    else:
        header_lineP = '\t'.join(str(e) for e in header_line)
        output.write('%s\n' % header_lineP)  

    for line in datafile:
        words = line.split()
        ch = int(words[0].split('chr')[-1])
        pos = int(round(float(words[1]),-1))
        content = words[2:]

        # read the input 2 as necessary
        while ((ch == inputFile2_CHR and pos > inputFile2_POS) or \
                (ch > inputFile2_CHR)) and (inputFile2_words != []):
                inputFile2_words = inputFile2.readline().split()
                try:
                    inputFile2_CHR = int(inputFile2_words[0].split('chr')[-1])
                    inputFile2_POS = int(round(float(inputFile2_words[1]), -1))
                    inputFile2_content = inputFile2_words[2:]
                except:
                    continue

        # find overlap
        if ch == inputFile2_CHR and pos == inputFile2_POS:
            content = inputFile2_content
            inputFile2_words = inputFile2.readline().split()
            if inputFile2_words != []:
                inputFile2_CHR = int(inputFile2_words[0].split('chr')[-1])
                inputFile2_POS = int(round(float(inputFile2_words[1]), -1))
                inputFile2_content = inputFile2_words[2:]          
  
        contentP =  '\t'.join(str(e) for e in content)
        output.write('chr%s\t%s\t%s\n' % (ch, pos, contentP))

        counter = calls.lineCounter(counter)

datafile.close()
inputFile2.close()
output.close()

print('Done!')
