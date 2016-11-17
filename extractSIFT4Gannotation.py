#! /usr/bin/env python

'''
This script extracts the SIFT4G annotation (http://sift.bii.a-star.edu.sg/sift4g/AnnotateVariants.html) for a given set of samples according to their genotypes.

http://sift.bii.a-star.edu.sg/sift4g/AnnotateVariants.html

# tab file:

CHROM	POS	REF	sample1	sample2	sample3	sample4
chr_1	2068	T	A	N	T	A
chr_1	2069	A	C	C	T	N

# SIFT4G annotation:

chr_1	2068	T	A	gene1	gene1.g	gene1.g	CDS	SYNONYMOUS	I	I	167	1.000	3.31	85	novel	TOLERATED
chr_1	2068	T	T	gene1	gene1.g	gene1.g	CDS	SYNONYMOUS	I	I	167	1.000	3.31	85	ref	TOLERATED
chr_1	2068	T	G	gene1	gene1.g	gene1.g	CDS	SYNONYMOUS	I	I	167	1.000	3.31	85	novel	TOLERATED
chr_1	2068	T	C	gene1	gene1.g	gene1.g	CDS	NONSYNONYMOUS	I	M	167	0.001	3.31	85	novel	DELETERIOUS
chr_1	2069	A	A	gene1	gene1.g	gene1.g	CDS	SYNONYMOUS	I	I	167	1.000	3.31	85	ref	TOLERATED
chr_1	2069	A	T	gene1	gene1.g	gene1.g	CDS	NONSYNONYMOUS	I	K	167	0.000	3.31	85	novel	DELETERIOUS
chr_1	2069	A	G	gene1	gene1.g	gene1.g	CDS	NONSYNONYMOUS	I	T	167	0.001	3.31	85	novel	DELETERIOUS
chr_1	2069	A	C	gene1	gene1.g	gene1.g	CDS	NONSYNONYMOUS	I	R	167	0.000	3.31	85	novel	DELETERIOUS

# output:

CHROM   POS sample1 sample2 sample3 sample4
chr_1   2068    1.000|TOLERATED NA  1.000|TOLERATED 1.000|TOLERATED
chr_1   2069    0.000|DELETERIOUS   0.000|DELETERIOUS   0.000|DELETERIOUS   NA

Note! Chromosome number must be separated by _.
For example, chr_1 - correct, chr1 - incorrect.

# command:

$ python extractSIFT4Gannotation.py -i test.sift4g -t test.tab -o test.output -s "sample1,sample2,sample3,sample4" -f "SIFT_SCORE,PREDICTION"

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-t', '--tab', help = 'tab delimited genotype file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-f', '--fields', help = 'annotation fields to extract. Possible options: CHROM, POS, REF_ALLELE, ALT_ALLELE, TRANSCRIPT_ID, GENE_ID, GENE_NAME, REGION, VARIANT_TYPE, REF_AMINO, ALT_AMINO, AMINO_POS, SIFT_SCORE, SIFT_MEDIAN, NUM_SEQS, dbSNP, SIFT_PREDICTION', type=str, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process (optional)', type=str, required=False)
args = parser.parse_args()

# check if samples names are given and if all sample names are present in a header
sampleNames = calls.checkSampleNames(args.samples, args.tab)

############################# program #############################

print('Opening the file...')

counter = 0

siftFile = open(args.input, 'r')
annotOptions = siftFile.readline().split()
fieldsIndex = calls.indexSamples([args.fields], annotOptions)

sift_words = siftFile.readline().split()
sift_chr = int(sift_words[0].split('_')[1])
sift_pos = int(sift_words[1])

with open(args.tab) as datafile:
  header_words = datafile.readline().split()

  # index samples
  sampCol = calls.indexSamples(sampleNames, header_words)

  # make output header
  print('Creating the output file...')
  output = open(args.output, 'w')
  ouput_header = header_words[0:2] + calls.selectSamples(sampCol, header_words)
  ouput_headerP = '\t'.join(str(el) for el in ouput_header)
  output.write('%s\n' % ouput_headerP)

############################### perform counting ####################

  for line in datafile:
    words = line.split()
    ch = int(words[0].split('_')[1])
    pos = int(words[1])

    # select samples
    tab_charaters = calls.selectSamples(sampCol, words)

    # find overlap in genomic position
    while (ch > sift_chr) or (ch == sift_chr and pos > sift_pos):
      sift_words = siftFile.readline().split()
      if sift_words == []:
        break
      sift_chr = int(sift_words[0].split('_')[1])
      sift_pos = int(sift_words[1])

    # extract all annotations of the position
    sift_alt = []
    sift_annot = []
    while (ch == sift_chr and pos == sift_pos):
      sift_alt.append(sift_words[3])
      sift_annot.append([sift_words[a] for a in fieldsIndex])
      #print words[0:2], sift_words[0:4] # for debugging
      sift_words = siftFile.readline().split()
      if sift_words == []:
        break
      sift_chr = int(sift_words[0].split('_')[1])
      sift_pos = int(sift_words[1])

    # find overlap in tab genotypes and sift ALT
    tab_annot = []
    for s in tab_charaters:
      if s in sift_alt: # check if genotype is A,T,G,C (not N)
        annotIndex = sift_alt.index(s)
        s_annot = sift_annot[annotIndex]
        s_annotP = '|'.join(str(el) for el in s_annot)
        tab_annot.append(s_annotP)
      else:
        tab_annot.append('NA')

    # output annotation
    if any(x != 'NA' for x in tab_annot):
      tab_annotP = '\t'.join(str(el) for el in tab_annot)
      chr_posP = '\t'.join(str(el) for el in words[0:2])
      output.write('%s\t%s\n' % (chr_posP, tab_annotP))

    # track progress
    counter += 1
    if counter % 100000 == 0:
      print str(counter), "lines processed"

datafile.close()
siftFile.close()
output.close()

print('Done!')
