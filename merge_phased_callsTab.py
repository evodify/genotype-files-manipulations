#!/usr/bin/python2

"""
This script merges phased sites into two-character coded genotype file.

# input:

CHROM   POS sample2_1   sample2_2   sample7_1   sample7_2 
chr_1   1   N   N   N   N
chr_1   2   T   C   C   T
chr_1   3   -   -   C   C
chr_1   4   T   T   T   T
chr_1   6   C   C   C   C
chr_2   1   A   A   A   A
chr_2   2   C   C   C   C
chr_2   3   N   N   N   N
chr_2   4   T   T   C   C
chr_2   5   C   C   T   T
chr_3   1   N   N   N   N
chr_3   2   T   A   C   C
chr_3   3   N   N   N   N
chr_3   4   T   T   T   T
chr_3   5   N   N   C   C
chr_4   1   N   N   C   C
chr_4   2   N   N   G   C

# output:

CHROM   POS sample2 sample7
chr_1   1   N/N N/N
chr_1   2   T/C C/T
chr_1   3   -/- C/C
chr_1   4   T/T T/T
chr_1   6   C/C C/C
chr_2   1   A/A A/A
chr_2   2   C/C C/C
chr_2   3   N/N N/N
chr_2   4   T/T C/C
chr_2   5   C/C T/T
chr_3   1   N/N N/N
chr_3   2   T/A C/C
chr_3   3   N/N N/N
chr_3   4   T/T T/T
chr_3   5   N/N C/C
chr_4   1   N/N C/C
chr_4   2   N/N G/C


# command:

$ python merge_phased_callsTab.py -i inputfile -o outputfile -s "sample2_1,sample2_2,sample7_1,sample7_2"

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process. Order matters. Look inside the script for details ', type=str, required=False)
args = parser.parse_args()

# check if samples names are given and if all sample names are present in a header
if args.samples:
  sampleNames = calls.checkSampleNames(args.samples, args.input)
else:
  raise IOError('Spesify samples names with the option -s')

############################# program #############################


counter = 0

print('Opening the file...\n')

with open(args.input) as datafile:
  header_words = datafile.readline().split()

  ChrPos = header_words[0:2]
  ChrPosP = '\t'.join(str(e) for e in ChrPos)

  # index samples
  sampCol = calls.indexSamples(sampleNames, header_words)

  # create lists for output
  sampColnames = calls.selectSamples(sampCol, header_words)

  # create merged column names
  idsheader = []
  for i in range(len(sampColnames)):
    if i % 2 == 0:
      name1 = sampColnames[i].split("_")[0]
      name2 = sampColnames[i+1].split("_")[0]
      if name1 == name2:
        idsheader.append(name1)
      else:
        raise KeyError("Sample is not paired. Sample name %s doesn't equal to sample name %s" % (name1,name2))
      header = '\t'.join(str(e) for e in idsheader)

  print('Creating the output file...\n')
  outname = '%s' % args.output
  outfile = open(outname, 'w')
  outfile.write('%s\t%s\n' % (ChrPosP, header))  

  print('Merging alleles...\n')

  for line in datafile:
    words = line.split()
    chr_pos = words[0:2]

    # create lists for output
    alleles = calls.selectSamples(sampCol, words)

    # check if all genotypes are correct
    calls.if_all_gt_correct(alleles, line)

    # merge
    mergedAlles = []
    for i in range(len(alleles)):
      if i % 2 == 0:
        gt = alleles[i]+"/"+alleles[i+1]
        mergedAlles.append(gt)
 
    chromPosP = '\t'.join(str(e) for e in chr_pos)
    mergedAllesP = '\t'.join(str(e) for e in mergedAlles)
    outfile.write("%s\t%s\n" % (chromPosP, mergedAllesP))
    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
outfile.close()
print('\nDone!')
