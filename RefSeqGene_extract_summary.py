#!/usr/bin/python2

"""
This script extracts summary info from NCBI Gene annotation available at ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/

The code was modified from https://www.biostars.org/p/2144/

# input refseqgene.genomic.gbff:

//
LOCUS       NG_027818             216459 bp    DNA     linear   PRI 31-MAR-2017
DEFINITION  Homo sapiens myeloid/lymphoid or mixed-lineage leukemia;
            translocated to, 10 (MLLT10), RefSeqGene on chromosome 10.
ACCESSION   NG_027818
VERSION     NG_027818.1
KEYWORDS    RefSeq; RefSeqGene.
SOURCE      Homo sapiens (human)
  ORGANISM  Homo sapiens
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
            Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini;
            Catarrhini; Hominidae; Homo.
COMMENT     REVIEWED REFSEQ: This record has been curated by NCBI staff. The
            reference sequence was derived from AL358780.22, AL161799.19,
            AL357372.12 and AL359697.23.
            This sequence is a reference standard in the RefSeqGene project.

            Summary: This gene encodes a transcription factor and has been
            identified as a partner gene involved in several chromosomal
            rearrangements resulting in various leukemias. Multiple transcript
            variants encoding different isoforms have been found for this gene.
            [provided by RefSeq, Sep 2010].
PRIMARY     REFSEQ_SPAN         PRIMARY_IDENTIFIER PRIMARY_SPAN        COMP
            1-51638             AL358780.22        60674-112311
            51639-167914        AL161799.19        23721-139996
            167915-209467       AL357372.12        2001-43553
            209468-216459       AL359697.23        2001-8992
FEATURES             Location/Qualifiers
     source          1..216459
                     /organism="Homo sapiens"
                     /mol_type="genomic DNA"
                     /db_xref="taxon:9606"
                     /chromosome="10"
                     /map="10p12.31"
     variation       5
                     /replace="c"
                     /replace="t"
                     /db_xref="dbSNP:1004761028"

# input gene_RefSeqGene:

#tax_id	GeneID	Symbol	RSG
9606	23209	MLC1	NG_009162.1
9606	4291	MLF1	NG_027720.1
9606	4292	MLH1	NG_007109.2
9606	27030	MLH3	NG_008649.1
9606	8028	MLLT10	NG_027818.1
9606	2862	MLNR	NG_029790.1
9606	79083	MLPH	NG_007286.1
9606	6945	MLX	NG_029442.1

# output :


# command:

$ python RefSeqGene_extract_summary.py -i inputfile -g geneID.refSeq -o outputfile

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""
############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the gbff input file', type=str, required=True)
parser.add_argument('-g', '--geneID', help = 'name of the gene_RefSeqGene file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
args = parser.parse_args()

############################# program #############################


# create gene ID dict
idDict = {}
with open(args.geneID) as idfile:
  ideHeader = idfile.readline()
  for idline in idfile:
    idwords = idline.split("\t")
    RSG = idwords[3].split(".")[0]
    if RSG in idDict.keys():
      idDict[RSG].append(idwords[2])
    else:
      idDict[RSG] = [idwords[2]]
idfile.close()

fout = open(args.output, 'w')
with open(args.input) as f:
  locus2comment = {}
  in_comment=False
  for line in f:
    if line[0:5] == "LOCUS":
      locus = line.split()[1]
      comment = ""
    elif line[0:7] == "COMMENT":
      in_comment=True
      comment += line.split("    ")[1].replace("\n", " ")
    elif line[0:7] == "PRIMARY":
      in_comment = False
      try:
        locus2comment[locus] = comment.split("Summary:")[1]
      except:
        locus2comment[locus] = comment
    elif in_comment:
      comment += line.split("            ")[1].replace("\n", " ")
  for locus in sorted(locus2comment):
    genesP = ','.join(str(e) for e in idDict[locus])
    fout.write(genesP + '\t' + locus2comment[locus] + '\n')

fout.close()
f.close()
print('Done!')
