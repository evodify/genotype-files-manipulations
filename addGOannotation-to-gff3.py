#! /usr/bin/env python2

'''
This script adds GO annotation to the gff3 file.

#Example input (gff3 or any other format that included the attribute ID=gene:):

1	247829	322180	ID=gene:ENSCAFG00000000001	Name=ENPP1
1	270005	275409	ID=gene:ENSCAFG00000030108
1	363607	364548	ID=gene:ENSCAFG00000000002
1	509125	565905	ID=gene:ENSCAFG00000000005	Name=PARD6G

#Example of GO annotation file

ENSCAFG00000000001	GO:0001889, GO:0003676, GO:0003824, GO:0004527, GO:0004528, GO:0004528, GO:0004551, GO:0004551, GO:0005044, GO:0005158, GO:0005158, GO:0005509, GO:0005515, GO:0005524, GO:0005524, GO:0005615, GO:0005615, GO:0005886, GO:0005886, GO:0005887, GO:0006091, GO:0006091, GO:0006796, GO:0006898, GO:0006955, GO:0008270, GO:0009143, GO:0009143, GO:0009612, GO:0009986, GO:0009986, GO:0016787, GO:0021549, GO:0021756, GO:0021766, GO:0021772, GO:0021987, GO:0030247, GO:0030308, GO:0030308, GO:0030505, GO:0030505, GO:0030643, GO:0030643, GO:0030730, GO:0030730, GO:0031953, GO:0031953, GO:0032869, GO:0032869, GO:0042803, GO:0042803, GO:0042995, GO:0043005, GO:0043025, GO:0045202, GO:0045599, GO:0045599, GO:0045719, GO:0045719, GO:0046034, GO:0046325, GO:0046325, GO:0046627, GO:0046627, GO:0046849, GO:0046872, GO:0047429, GO:0047429, GO:0050427, GO:0050427, GO:0071260, GO:0071320, GO:0071468, GO:0071560, GO:0090305, GO:0097440, GO:0097441
ENSCAFG00000000005	GO:0005080, GO:0005634, GO:0005938, GO:0007098, GO:0007163, GO:0007163, GO:0016324, GO:0017048, GO:0060341
ENSCAFG00000000007	GO:0003676, GO:0003677, GO:0030182, GO:0030307, GO:0034599, GO:0060548, GO:0071300
ENSCAFG00000024219	GO:0006364
ENSCAFG00000000008	GO:0000398, GO:0000398, GO:0005634, GO:0005682, GO:0005829, GO:0031965, GO:0046540, GO:0046540, GO:0046540, GO:0071005

#Example output:

1	247829	322180	ID=gene:ENSCAFG00000000001	Name=ENPP1	GO:0001889, GO:0003676, GO:0003824, GO:0004527, GO:0004528, GO:0004528, GO:0004551, GO:0004551, GO:0005044, GO:0005158, GO:0005158, GO:0005509, GO:0005515, GO:0005524, GO:0005524, GO:0005615, GO:0005615, GO:0005886, GO:0005886, GO:0005887, GO:0006091, GO:0006091, GO:0006796, GO:0006898, GO:0006955, GO:0008270, GO:0009143, GO:0009143, GO:0009612, GO:0009986, GO:0009986, GO:0016787, GO:0021549, GO:0021756, GO:0021766, GO:0021772, GO:0021987, GO:0030247, GO:0030308, GO:0030308, GO:0030505, GO:0030505, GO:0030643, GO:0030643, GO:0030730, GO:0030730, GO:0031953, GO:0031953, GO:0032869, GO:0032869, GO:0042803, GO:0042803, GO:0042995, GO:0043005, GO:0043025, GO:0045202, GO:0045599, GO:0045599, GO:0045719, GO:0045719, GO:0046034, GO:0046325, GO:0046325, GO:0046627, GO:0046627, GO:0046849, GO:0046872, GO:0047429, GO:0047429, GO:0050427, GO:0050427, GO:0071260, GO:0071320, GO:0071468, GO:0071560, GO:0090305, GO:0097440, GO:0097441
1	270005	275409	ID=gene:ENSCAFG00000030108	no GO annotation
1	363607	364548	ID=gene:ENSCAFG00000000002	no GO annotation
1	509125	565905	ID=gene:ENSCAFG00000000005	Name=PARD6G	GO:0005080, GO:0005634, GO:0005938, GO:0007098, GO:0007163, GO:0007163, GO:0016324, GO:0017048, GO:0060341


#command:

$  python addGOannotation-to-gff3.py -i test.gff3 -o test.GO.gff3 -g GOannot.csv

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-g', '--GOannotation', help = 'name of the GO annotation file (see the script header for examples)',
                    type=str, required=True)

args = parser.parse_args()

############################# program #############################

output = open(args.output, 'w')

# prepare GO database
with open(args.GOannotation) as GOfile:
  GOwords = GOfile.readline().split("\t")
  GOdicts = {GOwords[0]: GOwords[1].rstrip()}
  for GOline in GOfile:
    GOwords = GOline.split("\t")
    GOdicts[GOwords[0]] = GOwords[1].rstrip()
GOfile.close()

with open(args.input) as datafile:
  for line in datafile:
    if not line.startswith('#'): # header filter
      words = line.split("\t")
      for w in words:
        if 'ID=gene:' in w:
          geneName = w.replace("ID=gene:", "").rstrip()
      try:
        GOextract = GOdicts[geneName]
      except:
        GOextract = "NA"
      wordsP = '\t'.join(str(e) for e in words)
      output.write('%s\t%s\n' % (wordsP.rstrip(), GOextract))
      columns = []

datafile.close()
output.close()

print('Done!')
