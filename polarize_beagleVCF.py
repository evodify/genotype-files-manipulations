#! /usr/bin/env python

"""
This script polarizes BEAGLE phased genotypes relative to an outgroup/ancestral sequence.

After polarization:

0 - ancestral
1 - derived
The REF and ALT columns are not changed, but the ancestral allele is added in the filter field.

# input:

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  AM001   AM002
chr6    215     .       T       G       .       PASS    .       GT      1|0     1|0
chr6    245     .       A       C       .       PASS    .       GT      0|0     0|0
chr6    261     .       C       A       .       PASS    .       GT      0|0     0|0
chr6    280     .       C       A       .       PASS    .       GT      0|0     1|0
chr6    321     .       G       C       .       PASS    .       GT      0|0     0|0
chr6    328     .       C       A       .       PASS    .       GT      0|0     0|0
chr6    345     .       G       A       .       PASS    .       GT      0|0     0|0
chr6    365     .       A       G       .       PASS    .       GT      0|0     0|0
chr6    367     .       G       C       .       PASS    .       GT      0|0     0|0


# ancestral:

CHROM   POS     REF
chr6    215     T
chr6    245     G
chr6    261     N
chr6    280     C
chr6    321     G
chr6    345     A
chr6    367     G

# output:

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	AM001	AM002
chr6	215	.	T	G	.	ancestral=T	.	GT	1|0	1|0
chr6	245	.	A	C	.	ancestral=C	.	GT	1|1	1|1
chr6	280	.	C	A	.	ancestral=C	.	GT	0|0	1|0
chr6	321	.	G	C	.	ancestral=G	.	GT	0|0	0|0
chr6	345	.	G	A	.	ancestral=A	.	GT	1|1	1|1
chr6	367	.	G	C	.	ancestral=G	.	GT	0|0	0|0

# command:

$  python polarize_beagleVCF.py -i beagle.vcf -a ancestor.tab -o beagle.polarized

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""
############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the reference input file', type=str, required=True)
parser.add_argument('-a', '--ancestral', help = 'name of the outgroup/ancestral sequence file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-c', '--skip_checking_if_input_sorted', help = 'skip checking if input files are naturally sorted in chr and pos columns (type anything)', type=bool, required=False)

args = parser.parse_args()

if not args.skip_checking_if_input_sorted:
  calls.check_if_chr_pos_sorted(args.input)
  calls.check_if_chr_pos_sorted(args.ancestral)

############################# program #############################

counter = 0

# read the header
ances = open(args.ancestral, 'r')
ances_words = ances.readline()

output = open(args.output, 'w')

print('Opening the file...')
with open(args.input) as datafile:
  # make a header of the output file
  header_words = datafile.readline().split()
  print('Creating the output file...')
  header_wordsP = '\t'.join(str(e) for e in header_words)
  output.write("%s\n" % header_wordsP)

  # read the second line of the ancestral file
  ances_words = ances.readline().split()
  ances_ch = ances_words[0]
  ances_pos = int(ances_words[1])
  ances_gt = ances_words[2]

  for line in datafile:
    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

    words = line.split()
    ch = words[0]
    pos = int(words[1])
    samples_gt = words[9:]
    ref_gt = words[3]
    alt_gt = words[4]
    gt_info = words[0:9]

    # find overlap
    while (ch != ances_ch) or (ch == ances_ch and pos > ances_pos):
      ances_words = ances.readline().split()
      if ances_words == []:
        break
      else:
        ances_ch = ances_words[0]
        ances_pos = int(ances_words[1])
        ances_gt = ances_words[2]

    # introduce missing data if there is no overlap
    if pos != ances_pos:
      continue  # skip all missing data lines
    elif ances_gt == 'N':
      continue  # skip missing ancestral positions
    else: # polarize
      #gt_info[6] = "ancestral=" + ances_gt + "," + samples_gt[0] +","+ samples_gt[1]
      for i in range(len(samples_gt)):
        gt_hapl = samples_gt[i].split('|')
        if ref_gt != ances_gt and alt_gt == ances_gt:
          for j in range(len(gt_hapl)):
            if gt_hapl[j] == '0':
              gt_hapl[j] = '1'
            elif gt_hapl[j] == '1':
              gt_hapl[j] = '0'
        samples_gt[i] = '|'.join(str(e) for e in gt_hapl)

      gt_info[6] = "ancestral="+ ances_gt
      gt_infoP = '\t'.join(str(e) for e in gt_info)
      samples_gtP = '\t'.join(str(e) for e in samples_gt)
      output.write('%s\t%s\n' % (gt_infoP, samples_gtP))

datafile.close()
ances.close()
output.close()

print('Done!')
