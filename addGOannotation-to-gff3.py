#! /usr/bin/env python2

'''
This script extracts various info from the gff3 file.

#Example input:
1	Ensembl	chromosome	1	122678785	.	.	.	ID=chromosome:1;Alias=CM000001.3,NC_006583.3
1	tRNAscan	biological_region	212486	212568	21.5	+	.	external_name=Pseudo;logic_name=trnascan
1	ensembl	gene	247829	322180	.	-	.	ID=gene:ENSCAFG00000000001;Name=ENPP1;biotype=protein_coding;description=ectonucleotide pyrophosphatase/phosphodiesterase 1 [Source:VGNC Symbol%3BAcc:VGNC:40374];gene_id=ENSCAFG00000000001;logic_name=ensembl;version=3
1	ensembl	mRNA	247829	322180	.	-	.	ID=transcript:ENSCAFT00000000001;Parent=gene:ENSCAFG00000000001;Name=ENPP1-201;biotype=protein_coding;transcript_id=ENSCAFT00000000001;version=3
1	ensembl	three_prime_UTR	247829	252391	.	-	.	Parent=transcript:ENSCAFT00000000001
1	ensembl	exon	247829	252562	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000325317;constitutive=1;ensembl_end_phase=-1;ensembl_phase=0;exon_id=ENSCAFE00000325317;rank=25;version=1
1	ensembl	CDS	252392	252562	.	-	0	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	254466	254628	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000025;constitutive=1;ensembl_end_phase=0;ensembl_phase=2;exon_id=ENSCAFE00000000025;rank=24;version=1
1	ensembl	CDS	254466	254628	.	-	1	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	256197	256329	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000024;constitutive=1;ensembl_end_phase=2;ensembl_phase=1;exon_id=ENSCAFE00000000024;rank=23;version=1
1	ensembl	CDS	256197	256329	.	-	2	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	257162	257242	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000023;constitutive=1;ensembl_end_phase=1;ensembl_phase=1;exon_id=ENSCAFE00000000023;rank=22;version=1
1	ensembl	CDS	257162	257242	.	-	2	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	261834	261963	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000022;constitutive=1;ensembl_end_phase=1;ensembl_phase=0;exon_id=ENSCAFE00000000022;rank=21;version=1
1	ensembl	CDS	261834	261963	.	-	0	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	263524	263678	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000020;constitutive=1;ensembl_end_phase=0;ensembl_phase=1;exon_id=ENSCAFE00000000020;rank=20;version=3
1	ensembl	CDS	263524	263678	.	-	2	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	265578	265629	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000312619;constitutive=1;ensembl_end_phase=1;ensembl_phase=0;exon_id=ENSCAFE00000312619;rank=19;version=1
1	ensembl	CDS	265578	265629	.	-	0	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	266459	266628	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000018;constitutive=1;ensembl_end_phase=0;ensembl_phase=1;exon_id=ENSCAFE00000000018;rank=18;version=3
1	ensembl	CDS	266459	266628	.	-	2	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	268251	268338	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000017;constitutive=1;ensembl_end_phase=1;ensembl_phase=0;exon_id=ENSCAFE00000000017;rank=17;version=1
1	ensembl	CDS	268251	268338	.	-	0	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	269066	269135	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000016;constitutive=1;ensembl_end_phase=0;ensembl_phase=2;exon_id=ENSCAFE00000000016;rank=16;version=1
1	ensembl	CDS	269066	269135	.	-	1	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	270484	270611	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000015;constitutive=1;ensembl_end_phase=2;ensembl_phase=0;exon_id=ENSCAFE00000000015;rank=15;version=1
1	ensembl	CDS	270484	270611	.	-	0	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	271765	271796	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000014;constitutive=1;ensembl_end_phase=0;ensembl_phase=1;exon_id=ENSCAFE00000000014;rank=14;version=1
1	ensembl	CDS	271765	271796	.	-	2	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	273565	273696	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000013;constitutive=1;ensembl_end_phase=1;ensembl_phase=1;exon_id=ENSCAFE00000000013;rank=13;version=1
1	ensembl	CDS	273565	273696	.	-	2	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	275918	276026	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000012;constitutive=1;ensembl_end_phase=1;ensembl_phase=0;exon_id=ENSCAFE00000000012;rank=12;version=1
1	ensembl	CDS	275918	276026	.	-	0	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	278245	278317	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000011;constitutive=1;ensembl_end_phase=0;ensembl_phase=2;exon_id=ENSCAFE00000000011;rank=11;version=1
1	ensembl	CDS	278245	278317	.	-	1	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	278628	278693	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000010;constitutive=1;ensembl_end_phase=2;ensembl_phase=2;exon_id=ENSCAFE00000000010;rank=10;version=1
1	ensembl	CDS	278628	278693	.	-	1	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	279595	279704	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000009;constitutive=1;ensembl_end_phase=2;ensembl_phase=0;exon_id=ENSCAFE00000000009;rank=9;version=1
1	ensembl	CDS	279595	279704	.	-	0	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	280770	280889	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000008;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=ENSCAFE00000000008;rank=8;version=1
1	ensembl	CDS	280770	280889	.	-	0	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	281894	281973	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000007;constitutive=1;ensembl_end_phase=0;ensembl_phase=1;exon_id=ENSCAFE00000000007;rank=7;version=1
1	ensembl	CDS	281894	281973	.	-	2	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	283862	283959	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000006;constitutive=1;ensembl_end_phase=1;ensembl_phase=2;exon_id=ENSCAFE00000000006;rank=6;version=1
1	ensembl	CDS	283862	283959	.	-	1	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	285904	285964	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000005;constitutive=1;ensembl_end_phase=2;ensembl_phase=1;exon_id=ENSCAFE00000000005;rank=5;version=2
1	ensembl	CDS	285904	285964	.	-	2	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	286842	286967	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000004;constitutive=1;ensembl_end_phase=1;ensembl_phase=1;exon_id=ENSCAFE00000000004;rank=4;version=2
1	ensembl	CDS	286842	286967	.	-	2	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	287787	287903	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000003;constitutive=1;ensembl_end_phase=1;ensembl_phase=1;exon_id=ENSCAFE00000000003;rank=3;version=1
1	ensembl	CDS	287787	287903	.	-	2	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	289674	289746	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000000002;constitutive=1;ensembl_end_phase=1;ensembl_phase=0;exon_id=ENSCAFE00000000002;rank=2;version=1
1	ensembl	CDS	289674	289746	.	-	0	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	CDS	321851	322063	.	-	0	ID=CDS:ENSCAFP00000000001;Parent=transcript:ENSCAFT00000000001;protein_id=ENSCAFP00000000001
1	ensembl	exon	321851	322180	.	-	.	Parent=transcript:ENSCAFT00000000001;Name=ENSCAFE00000236366;constitutive=1;ensembl_end_phase=0;ensembl_phase=-1;exon_id=ENSCAFE00000236366;rank=1;version=2
1	ensembl	five_prime_UTR	322064	322180	.	-	.	Parent=transcript:ENSCAFT00000000001
1	tRNAscan	biological_region	254125	254207	22.4	-	.	external_name=Pseudo;logic_name=trnascan
1	ensembl	gene	270005	275409	.	+	.	ID=gene:ENSCAFG00000030108;biotype=protein_coding;gene_id=ENSCAFG00000030108;logic_name=ensembl;version=1
1	ensembl	mRNA	270005	275409	.	+	.	ID=transcript:ENSCAFT00000043967;Parent=gene:ENSCAFG00000030108;biotype=protein_coding;transcript_id=ENSCAFT00000043967;version=1
1	ensembl	exon	270005	270109	.	+	.	Parent=transcript:ENSCAFT00000043967;Name=ENSCAFE00000311607;constitutive=1;ensembl_end_phase=0;ensembl_phase=-1;exon_id=ENSCAFE00000311607;rank=1;version=1
1	ensembl	CDS	270005	270109	.	+	0	ID=CDS:ENSCAFP00000042539;Parent=transcript:ENSCAFT00000043967;protein_id=ENSCAFP00000042539
1	ensembl	CDS	274549	275367	.	+	0	ID=CDS:ENSCAFP00000042539;Parent=transcript:ENSCAFT00000043967;protein_id=ENSCAFP00000042539
1	ensembl	exon	274549	275409	.	+	.	Parent=transcript:ENSCAFT00000043967;Name=ENSCAFE00000285925;constitutive=1;ensembl_end_phase=-1;ensembl_phase=0;exon_id=ENSCAFE00000285925;rank=2;version=1
1	ensembl	three_prime_UTR	275368	275409	.	+	.	Parent=transcript:ENSCAFT00000043967
1	tRNAscan	biological_region	305756	305837	22.5	+	.	external_name=Pseudo;logic_name=trnascan
1	cpg	biological_region	321600	322578	1.25e+03	.	.	external_name=oe %3D 0.83;logic_name=cpg
1	Eponine	biological_region	321651	321655	0.999	+	.	logic_name=eponine
1	.	biological_region	321695	322000	1	+	.	external_name=rank %3D 2;logic_name=firstef
1	.	biological_region	321696	322585	1	+	.	external_name=rank %3D 1;logic_name=firstef
1	Eponine	biological_region	321793	321796	0.999	+	.	logic_name=eponine
1	.	biological_region	321851	322276	1	-	.	external_name=rank %3D 1;logic_name=firstef
1	Eponine	biological_region	321880	321883	0.999	-	.	logic_name=eponine
1	Eponine	biological_region	321935	321938	0.999	+	.	logic_name=eponine
1	Eponine	biological_region	321971	321973	0.999	+	.	logic_name=eponine
1	Eponine	biological_region	322036	322036	0.999	-	.	logic_name=eponine
1	Eponine	biological_region	322051	322054	0.999	-	.	logic_name=eponine
1	Eponine	biological_region	322074	322085	0.999	-	.	logic_name=eponine
1	Eponine	biological_region	322090	322093	0.999	+	.	logic_name=eponine
1	Eponine	biological_region	322130	322139	0.999	+	.	logic_name=eponine
1	Eponine	biological_region	322162	322167	0.999	-	.	logic_name=eponine
1	Eponine	biological_region	322222	322229	0.999	-	.	logic_name=eponine
1	Eponine	biological_region	322264	322269	0.999	-	.	logic_name=eponine
1	Eponine	biological_region	322517	322521	0.999	-	.	logic_name=eponine
1	cpg	biological_region	329739	330203	184	.	.	external_name=oe %3D 1.00;logic_name=cpg
1	ensembl	gene	363607	364548	.	+	.	ID=gene:ENSCAFG00000000002;biotype=protein_coding;gene_id=ENSCAFG00000000002;logic_name=ensembl;version=3
1	ensembl	mRNA	363607	364548	.	+	.	ID=transcript:ENSCAFT00000000003;Parent=gene:ENSCAFG00000000002;biotype=protein_coding;transcript_id=ENSCAFT00000000003;version=3
1	ensembl	exon	363607	364548	.	+	.	Parent=transcript:ENSCAFT00000000003;Name=ENSCAFE00000000028;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSCAFE00000000028;rank=1;version=3
1	ensembl	CDS	363607	364548	.	+	0	ID=CDS:ENSCAFP00000041865;Parent=transcript:ENSCAFT00000000003;protein_id=ENSCAFP00000041865
1	.	biological_region	418153	418238	0.809	+	.	external_name=rank %3D 1;logic_name=firstef
1	tRNAscan	biological_region	418909	418991	22.5	-	.	external_name=Pseudo;logic_name=trnascan
1	cpg	biological_region	423322	423803	257	.	.	external_name=oe %3D 0.73;logic_name=cpg
1	tRNAscan	biological_region	470353	470436	25.2	-	.	external_name=Pseudo;logic_name=trnascan
1	tRNAscan	biological_region	474050	474131	36.4	+	.	external_name=SeC;logic_name=trnascan
1	.	biological_region	487953	488099	1	+	.	external_name=rank %3D 1;logic_name=firstef
1	tRNAscan	biological_region	507349	507430	22.8	-	.	external_name=Pseudo;logic_name=trnascan
1	ensembl	gene	509125	565905	.	+	.	ID=gene:ENSCAFG00000000005;Name=PARD6G;biotype=protein_coding;description=par-6 family cell polarity regulator gamma [Source:HGNC Symbol%3BAcc:HGNC:16076];gene_id=ENSCAFG00000000005;logic_name=ensembl;version=3
1	ensembl	mRNA	509125	565905	.	+	.	ID=transcript:ENSCAFT00000000006;Parent=gene:ENSCAFG00000000005;Name=PARD6G-201;biotype=protein_coding;transcript_id=ENSCAFT00000000006;version=3
1	ensembl	exon	509125	509139	.	+	.	Parent=transcript:ENSCAFT00000000006;Name=ENSCAFE00000317723;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=ENSCAFE00000317723;rank=1;version=1
1	ensembl	CDS	509125	509139	.	+	0	ID=CDS:ENSCAFP00000000006;Parent=transcript:ENSCAFT00000000006;protein_id=ENSCAFP00000000006
1	ensembl	exon	509347	509403	.	+	.	Parent=transcript:ENSCAFT00000000006;Name=ENSCAFE00000228100;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=ENSCAFE00000228100;rank=2;version=2



#Example output:

1	247829	322180	ID=gene:ENSCAFG00000000001	Name=ENPP1
1	270005	275409	ID=gene:ENSCAFG00000030108
1	363607	364548	ID=gene:ENSCAFG00000000002
1	509125	565905	ID=gene:ENSCAFG00000000005	Name=PARD6G

#command:

$ python annotate_GO_withStatFromSlidingWindowsStats.py -i input.table -o output.tab -s Fst -t GOref.csv

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-c', '--columns', help = 'comma separated names of columns to output', type=str, required=True)
parser.add_argument('-t', '--type', help = 'type of feature e.g. gene, mRNA, exon, CDS', type=str, required=False)
parser.add_argument('-a', '--attributes', help = 'comma separated names of features to extract from annotation e.g. ID, Name, Alias.) ', type=str, required=False)

args = parser.parse_args()

############################# program #############################

output = open(args.output, 'w')

colExtract = args.columns.split(",")
for col in colExtract:
  if int(col) > 9:
    raise IOError('There is only 9 columns. The column %s will be empty' % col)

with open(args.input) as datafile:
  for line in datafile:
    if not line.startswith('#'): # header filter
      words = line.split("\t")
      if words[2] == args.type: # type filter
        columns = []
        for col in colExtract:
          columns.append(words[int(col)-1])
        attrWords = words[8].split(';')
        for att in args.attributes.split(','):
          for attrW in attrWords:
            if str(att)+'=' in str(attrW):
              columns.append(attrW)
        columnsP = '\t'.join(str(e) for e in columns)
        output.write('%s\n' % columnsP)
        columns = []

datafile.close()
output.close()

print('Done!')
