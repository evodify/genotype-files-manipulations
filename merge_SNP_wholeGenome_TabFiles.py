#! /usr/bin/env python

"""
This script merges whole genome and SNPs tab files. This is needed because non-polymorphic sites and SNPs are filtered differently with GATK.

Note! Chromosome number in the first column must be separated by _.
For example, chr_1 - correct, chr1 - incorrect.

# whole_Genome.tab:

CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   1   A   A   A   A   A   A   A   A   A
chr_1   2   N   C   C   T   C   C   T   C   C
chr_1   3   T   A   A   T   T   C   T   A   T
chr_1   4   C   C   C   C   C   C   C   C   C
chr_1   5   C   C   T   C   C   C   C   C   C

# SNPs.tab:

CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   2   N   T   C   T   C   C   T   T   C
chr_1   4   A   C   A   C   A   C   A   C   C


# output:

CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   1   A   A   A   A   A   A   A   A   A
chr_1   2   N   T   C   T   C   C   T   T   C
chr_1   3   N   N   N   N   N   N   N   N   N
chr_1   4   N   C   A   C   A   C   A   C   C
chr_1   5   C   C   C   C   C   C   C   C   C


# command:

$ python merge_SNP_wholeGenome_TabFiles.py -g whole_Genome.tab -s SNPs.tab -o output.tab

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""
############################# modules #############################

import calls # my custom module

############################# options #############################


parser = calls.MyParser()
parser.add_argument('-g', '--whole_genome_input', help = 'name of the whole genome input file', type=str, required=True)
parser.add_argument('-s', '--SNPs_input', help = 'name of the SNPs file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
args = parser.parse_args()

############################# program #############################

counter = 0
WG_polymorphic = 0
# read the header of the larger file
SNPsFile = open(args.SNPs_input, 'r')
SNPsFile_header = SNPsFile.readline().split()

output = open(args.output, 'w')

print('Opening the file...')
with open(args.whole_genome_input) as WG:
  WG_header_line = WG.readline().split()
  if WG_header_line != SNPsFile_header: # check the input samples order
    raise IOError('Headers are different in the two files')

  print('Creating the output file...')
  WG_headerP = '\t'.join(str(e) for e in WG_header_line)
  output.write('%s\n' % WG_headerP)

  # read the second line of the larger file
  SNPsFile_words = SNPsFile.readline().split()
  SNPsFile_ch = int(SNPsFile_words[0].split('_')[1])
  SNPsFile_pos = int(SNPsFile_words[1])
  SNPsFile_gt = SNPsFile_words[2:]

  for line in WG:
    WG_words = line.split()
    WG_ch = int(WG_words[0].split('_')[1])
    WG_pos = int(WG_words[1])
    WG_gt = WG_words[2:]

    # find overlap
    if WG_ch == SNPsFile_ch and WG_pos == SNPsFile_pos:
      SNPsFile_gtP = '\t'.join(str(e) for e in (WG_words[0:2]+SNPsFile_gt))
      output.write('%s\n' % SNPsFile_gtP)
      SNPsFile_words = SNPsFile.readline().split()
      if SNPsFile_words == []:
        continue
      else:
        SNPsFile_ch = int(SNPsFile_words[0].split('_')[1])
        SNPsFile_pos = int(SNPsFile_words[1])
        SNPsFile_gt = SNPsFile_words[2:]
    else:
      if calls.is_polymorphic(WG_gt): # replace polymorphic sites in WG with Ns
        WG_gtN = ["N"] * len(WG_gt)
        WG_polymorphic += 1
        WG_gtP = '\t'.join(str(e) for e in (WG_words[0:2]+WG_gtN))
      else:
        WG_gtP = '\t'.join(str(e) for e in (WG_words[0:2]+WG_gt))
      output.write('%s\n' % WG_gtP)

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"
print "False polymorphic sites in the whole genome = "+str(WG_polymorphic)
WG.close()
SNPsFile.close()
output.close()

print('Done!')
