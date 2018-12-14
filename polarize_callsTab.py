#! /usr/bin/env python
"""
This script polarizes the genotype data relative to
an outgroup/ancestral sequence.

# input:

CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7
chr_1   1   A   W   N   N   A   N   N   N
chr_1   2   C   Y   Y   N   C   C   T   C
chr_1   3   C   N   -   N   C   C   C   C
chr_1   4   T   T   T   N   T   T   T   T
chr_1   6   C   N   C   N   C   C   C   C
chr_2   1   A   A   A   N   A   A   A   A
chr_2   2   C   C   C   N   C   C   C   C
chr_2   3   C   N   N   N   N   N   N   N
chr_2   4   C   C   T   C   C   C   C   C
chr_2   5   T   T   C   T   C   T   Y   T
chr_3   1   G   G   N   N   G   N   N   N
chr_3   2   C   S   C   N   C   C   N   C
chr_3   3   N   N   N   N   N   N   N   N
chr_3   4   N   T   T   N   T   T   T   T


# ancestral:

CHROM   POS ancestor
chr_1   1   W
chr_1   2   C
chr_1   3   C
chr_1   4   G
chr_1   6   C
chr_2   1   A
chr_2   2   T
chr_2   3   C
chr_2   4   C
chr_2   5   T
chr_3   1   G
chr_3   2   S
chr_3   3   N
chr_3   4   C


# output:

CHROM	POS	ANC	DER	REF	sample1	sample2	sample3	sample4	sample5	sample6	sample7
chr_1	1	A	T	0/0	1/0	./.	./.	0/0	./.	./.	./.
chr_1	2	C	T	0/0	1/0	1/0	./.	0/0	0/0	1/1	0/0
chr_1	3	C	-	0/0	./.	1/1	./.	0/0	0/0	0/0	0/0
chr_1	4	G	T	1/1	1/1	1/1	./.	1/1	1/1	1/1	1/1
chr_1	6	C	T	0/0	./.	0/0	./.	0/0	0/0	0/0	0/0
chr_2	1	A	T	0/0	0/0	0/0	./.	0/0	0/0	0/0	0/0
chr_2	2	T	C	1/1	1/1	1/1	./.	1/1	1/1	1/1	1/1
chr_2	3	C	C	0/0	./.	./.	./.	./.	./.	./.	./.
chr_2	4	C	T	0/0	0/0	1/1	0/0	0/0	0/0	0/0	0/0
chr_2	5	T	C	0/0	0/0	1/1	0/0	1/1	0/0	0/1	0/0
chr_3	1	G	C	0/0	0/0	./.	./.	0/0	./.	./.	./.
chr_3	2	C	G	0/0	0/1	0/0	./.	0/0	0/0	./.	0/0
chr_3	4	C	T	./.	1/1	1/1	./.	1/1	1/1	1/1	1/1


# command:

$ python polarize_callsTab.py \
    -i input.tab \
    -a ancestral.tab \
    -s "sample1,sample2,sample3,sample4,sample5,sample6,sample7" \
    -o output

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""
############################# modules #############################

import calls  # my custom module

def splitAncestral(ancestral):
    if ancestral in 'RYMKSW':
        ancestral = calls.pseudoPhase(ancestral)
    return ancestral[0]

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument(
    '-i',
    '--input',
    help='name of the reference input file',
    type=str,
    required=True)
parser.add_argument(
    '-a',
    '--ancestral',
    help='name of the outgroup/ancestral sequence file',
    type=str,
    required=True)
parser.add_argument(
    '-o', '--output',
    help='name of the output file',
    type=str,
    required=True)
parser.add_argument(
    '-s',
    '--samples',
    help='column names of the samples to process (optional)',
    type=str,
    required=False)
args = parser.parse_args()

sampleNames = calls.checkSampleNames(args.samples, args.input)

############################# program #############################

counter = 0

# read the header
ances = open(args.ancestral, 'r')
ances_words = ances.readline()

output = open(args.output, 'w')

print('Opening the file...')
with open(args.input) as datafile:
    header_words = datafile.readline().split()

    sampCol = calls.indexSamples(sampleNames, header_words)

    print('Creating the output file...')
    samples_head = calls.selectSamples(sampCol, header_words)
    samples_headP = '\t'.join(str(e) for e in samples_head)
    output.write('%s\t%s\t%s\t%s\t%s\n' % ('CHROM', 'POS', 'ANC',
                           'DER', samples_headP))

    # read the second line of the ancestral file
    ances_words = ances.readline().split()
    if '_' in ances_words[0]:
        ances_ch = int(ances_words[0].split('_')[1])
    else:
        ances_ch = int(ances_words[0].split('chr')[1])
    ances_pos = int(ances_words[1])
    ances_gt = splitAncestral(ances_words[2])

    for line in datafile:
        words = line.split()
        if '_' in words[0]:
            ch = int(words[0].split('_')[1])
        else:
            ch = int(words[0].split('chr')[1])
        pos = int(words[1])

        # find overlap
        while (ch > ances_ch) or (ch == ances_ch and pos > ances_pos):
            ances_words = ances.readline().split()
            if ances_words == []:
                break
            else:
                if '_' in ances_words[0]:
                    ances_ch = int(ances_words[0].split('_')[1])
                else:
                    ances_ch = int(ances_words[0].split('chr')[1])
                ances_pos = int(ances_words[1])
                ances_gt = splitAncestral(ances_words[2])

        if pos != ances_pos:
            continue  # skip all missing data lines
        else:
            samples_gt = calls.selectSamples(sampCol, words)
            if all('/' in gt for gt in samples_gt):
                samples_gtPhased = samples_gt
            else:
                samples_gtPhased =  calls.pseudoPhase(samples_gt)
            samples_gtPlarized = []
            derived = 'N'
            for i in range(len(samples_gtPhased)):
                gtS = samples_gtPhased[i].split('/')
                if ances_gt == 'N' or gtS[0] == 'N' or gtS[1] == 'N':
                    gtS = ['.', '.']
                else:
                    if (gtS[0] in ances_gt) and (gtS[1] not in ances_gt):
                        derived = gtS[1]
                        gtS[0] = '0'
                        gtS[1] = '1'
                    elif (gtS[0] not in ances_gt) and (gtS[1] in ances_gt):
                        derived = gtS[0]
                        gtS[0] = '1'
                        gtS[1] = '0'
                    elif (gtS[0] not in ances_gt) and (gtS[1] not in ances_gt):
                        if gtS[0] == gtS[1]:
                            derived = gtS[0]
                        else:
                            derived = '/'.join(str(e) for e in gtS)
                        gtS[0] = '1'
                        gtS[1] = '1'
                    elif (gtS[0] in ances_gt) and (gtS[1] in ances_gt) and
                        gtS[0] != gtS[1]:
                        gtS = random.choice([['0', '1'], ['1', '0']])
                    else:
                        raise IOError('Unexpected case: Ancestral:%s, GT:%s '
                                      'in %s %s' %
                                      (ances_gt, gtS, words[0], words[1]))
                gtSP = '/'.join(str(e) for e in gtS)
                samples_gtPlarized.append(gtSP)
            if all('./.' in gt for gt in samples_gtPlarized) or derived == 'N':
                continue
            else:
                wordsP = '\t'.join(str(e) for e in samples_gtPlarized)
        output.write('%s\t%s\t%s\t%s\t%s\n' % 
                           (words[0], words[1], ances_words[2], derived, wordsP))
        counter = calls.lineCounter(counter)
datafile.close()
ances.close()
output.close()
print('Done!')
