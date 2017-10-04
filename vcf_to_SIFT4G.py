#! /usr/bin/env python2

'''
This script converts a VCF file to SIFT4G input.

#Example input:

#CHROM  POS ID  REF
scaffold_1  1   .   C
scaffold_1  2   .   CAA
scaffold_1  3   .   T
scaffold_1  4   .   A
scaffold_1  5   .   A
scaffold_1  6   .   A
scaffold_1  7   .   C
scaffold_1  8   .   C
scaffold_1  9   .   C

#Example output:

#CHROM  POS ID  REF ALT QUAL    FILTER  INFO
scaffold_1  1   .   C   A   .   .   .
scaffold_1  1   .   C   T   .   .   .
scaffold_1  1   .   C   G   .   .   .
scaffold_1  1   .   C   C   .   .   .
scaffold_1  3   .   T   A   .   .   .
scaffold_1  3   .   T   T   .   .   .
scaffold_1  3   .   T   G   .   .   .
scaffold_1  3   .   T   C   .   .   .
scaffold_1  4   .   A   A   .   .   .
scaffold_1  4   .   A   T   .   .   .
scaffold_1  4   .   A   G   .   .   .
scaffold_1  4   .   A   C   .   .   .
scaffold_1  5   .   A   A   .   .   .
scaffold_1  5   .   A   T   .   .   .
scaffold_1  5   .   A   G   .   .   .
scaffold_1  5   .   A   C   .   .   .
scaffold_1  6   .   A   A   .   .   .
scaffold_1  6   .   A   T   .   .   .
scaffold_1  6   .   A   G   .   .   .
scaffold_1  6   .   A   C   .   .   .
scaffold_1  7   .   C   A   .   .   .
scaffold_1  7   .   C   T   .   .   .
scaffold_1  7   .   C   G   .   .   .
scaffold_1  7   .   C   C   .   .   .
scaffold_1  8   .   C   A   .   .   .
scaffold_1  8   .   C   T   .   .   .
scaffold_1  8   .   C   G   .   .   .
scaffold_1  8   .   C   C   .   .   .
scaffold_1  9   .   C   A   .   .   .
scaffold_1  9   .   C   T   .   .   .
scaffold_1  9   .   C   G   .   .   .
scaffold_1  9   .   C   C   .   .   .

#command:

$ python vcf_to_SIFT4G.py -i input.vcf -o output.vcf

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
args = parser.parse_args()

############################# program #############################

counter = 0

print('Opening the file...')

with open(args.input) as datafile:
  header_line = datafile.readline()

  print('Creating the output...')

  outfile = open(args.output, 'w')
  outfile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n') 

  print('Converting...')

  for line in datafile:
    words = line.split()
    chr_pos_ref = words[0:4]
    ref = words[3]

    if len(ref) == 1:  # to skip insertions
      for nucl in ['A', 'T', 'G', 'C']:
        chr_pos_refP = '\t'.join(str(e) for e in chr_pos_ref)
        outfile.write("%s\t%s\t.\t.\t.\n" % (chr_pos_refP, nucl))

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
outfile.close()

print('Done!')
 
