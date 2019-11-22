#! /usr/bin/env python
"""
Merges CNV tab file with the provided list of bins.

# CNV tab:

CHROM	POS	END	sample.CN
chr6	45894001	46949000	2
chr6	46949001	46956000	1
chr6	46956001	47139000	2
chr6	47139001	47148000	1
chr6	47148001	48052000	2
chr6	48052001	48053000	5 

# Interval bins:

CHROM	POS	END
chr6	45894001	45895000
chr6	45895001	46072000
chr6	46072001	46089000
chr6	46089001	46091000
chr6	46091001	46097000
chr6	46097001	46466000
chr6	46466001	46467000
chr6	46467001	46617000
chr6	46617001	46648000
chr6	46648001	46649000
chr6	46649001	46792000
chr6	46792001	46800000
chr6	46800001	46890000
chr6	46890001	46939000
chr6	46939001	46949000
chr6	46949001	46951000
chr6	46951001	46955000
chr6	46955001	46956000
chr6	46956001	46957000
chr6	46957001	47003000
chr6	47003001	47004000
chr6	47004001	47046000
chr6	47046001	47047000
chr6	47047001	47053000
chr6	47053001	47107000
chr6	47107001	47133000
chr6	47133001	47139000
chr6	47139001	47140000
chr6	47140001	47147000
chr6	47147001	47148000
chr6	47148001	47678000
chr6	47678001	47680000
chr6	47680001	47740000
chr6	47740001	47756000
chr6	47756001	47821000
chr6	47821001	47822000
chr6	47822001	47823000
chr6	47823001	47826000
chr6	47826001	48052000
chr6	48052001	48053000
chr6	48053001	48153000
chr6	48153001	48202000
chr6	48202001	48416000

# output:

CHROM	POS	END	sample.CN
chr6	45894001	45895000	2
chr6	45895001	46072000	2
chr6	46072001	46089000	2
chr6	46089001	46091000	2
chr6	46091001	46097000	2
chr6	46097001	46466000	2
chr6	46466001	46467000	2
chr6	46467001	46617000	2
chr6	46617001	46648000	2
chr6	46648001	46649000	2
chr6	46649001	46792000	2
chr6	46792001	46800000	2
chr6	46800001	46890000	2
chr6	46890001	46939000	2
chr6	46939001	46949000	2
chr6	46949001	46951000	1
chr6	46951001	46955000	1
chr6	46955001	46956000	1
chr6	46956001	46957000	2
chr6	46957001	47003000	2
chr6	47003001	47004000	2
chr6	47004001	47046000	2
chr6	47046001	47047000	2
chr6	47047001	47053000	2
chr6	47053001	47107000	2
chr6	47107001	47133000	2
chr6	47133001	47139000	2
chr6	47139001	47140000	1
chr6	47140001	47147000	1
chr6	47147001	47148000	1
chr6	47148001	47678000	2
chr6	47678001	47680000	2
chr6	47680001	47740000	2
chr6	47740001	47756000	2
chr6	47756001	47821000	2
chr6	47821001	47822000	2
chr6	47822001	47823000	2
chr6	47823001	47826000	2
chr6	47826001	48052000	2
chr6	48052001	48053000	5
chr6	48053001	48153000	NA
chr6	48153001	48202000	NA
chr6	48202001	48416000	NA


# command:

$ python merge_CNVs_tabs.py \
    -r cnv_bins.bed \
    -i cnv_sample.tab \
    -o cnv_sample_merged.tab

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

tabFile = open(args.input_to_merge, 'r')
tab_header = tabFile.readline()
tab_words = tabFile.readline().split()
tab_CHR = tab_words[0]
tab_startPOS = int(tab_words[1])
tab_endPOS = int(tab_words[2])
tab_CN = tab_words[3]
    
output = open(args.output, 'w')
output.write(tab_header)

print('Opening the file...')
with open(args.reference_input) as refFile:
    header_line = refFile.readline()

    for line in refFile:
        words = line.split()
        ch = words[0]
        startPOS = int(words[1])
        endPOS = int(words[2])

        # read the tab file as necessary
        while ((ch == tab_CHR and startPOS > tab_endPOS) or \
                (ch != tab_CHR)) and (tab_words != []):
                tab_words = tabFile.readline().split()
                try:
                    tab_CHR = tab_words[0]
                    tab_startPOS = int(tab_words[1])
                    tab_endPOS = int(tab_words[2])
                    tab_CN = tab_words[3]
                except:
                    continue

        # find a match
        if ch == tab_CHR and startPOS >= tab_startPOS and endPOS <= tab_endPOS:
            if tab_words != []:
                output.write('%s\t%s\t%s\t%s\n' % (ch, startPOS, endPOS, tab_CN))
        else:
            output.write('%s\t%s\t%s\tNA\n' % (ch, startPOS, endPOS))

        counter = calls.lineCounter(counter)

refFile.close()
tabFile.close()
output.close()

print('Done!')
