#! /usr/bin/env python2

'''
This script converts the Ensembl.dat file to the GO reference file used in the topGO R program.

#Example input:
ID   1    standard; DNA; HTG; 122678785 BP.
XX
AC   chromosome:CanFam3.1:1:1:122678785:1
XX
SV   1.CanFam3.1
XX
DT   4-SEP-2018
XX
DE   Canis lupus familiaris chromosome 1 CanFam3.1 full sequence 1..122678785
DE   annotated by Ensembl
XX
KW   .
XX
OS   Canis lupus familiaris (dog)
OC   Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia;
OC   Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi;
OC   Sarcopterygii; Dipnotetrapodomorpha; Tetrapoda; Amniota; Mammalia; Theria;
OC   Eutheria; Boreoeutheria; Laurasiatheria; Carnivora; Caniformia; Canidae;
OC   Canis lupus.
XX
CC   This sequence was annotated by Ensembl(www.ensembl.org). Please visit the
CC   Ensembl or EnsemblGenomes web site, http://www.ensembl.org/ or
CC   http://www.ensemblgenomes.org/ for more information.
XX
CC   All feature locations are relative to the first (5') base of the sequence
CC   in this file.  The sequence presented is always the forward strand of the
CC   assembly. Features that lie outside of the sequence contained in this file
CC   have clonal location coordinates in the format: <clone
CC   accession>.<version>:<start>..<end>
XX
CC   The /gene indicates a unique id for a gene, /note="transcript_id=..." a
CC   unique id for a transcript, /protein_id a unique id for a peptide and
CC   note="exon_id=..." a unique id for an exon. These ids are maintained
CC   wherever possible between versions.
XX
CC   All the exons and transcripts in Ensembl are confirmed by similarity to
CC   either protein or cDNA sequences.
XX
FH   Key             Location/Qualifiers
FT   source          1..122678785
FT                   /organism="Canis lupus familiaris"
FT                   /db_xref="taxon:9615"
FT   gene            722179..735934
FT                   /gene=ENSCAFG00000000008.3
FT                   /locus_tag="TXNL4A"
FT                   /note="thioredoxin like 4A [Source:VGNC
FT                   Symbol;Acc:VGNC:48019]"
FT   mRNA            join(722179..722324,722691..722877,731542..731645,
FT                   734838..735934)
FT                   /gene="ENSCAFG00000000008.3"
FT                   /standard_name="ENSCAFT00000000009.3"
FT   CDS             join(722725..722877,731542..731645,734838..735009)
FT                   /gene="ENSCAFG00000000008.3"
FT                   /protein_id="ENSCAFP00000000008.3"
FT                   /note="transcript_id=ENSCAFT00000000009.3"
FT                   /db_xref="RefSeq_mRNA_predicted:XM_005615276"
FT                   /db_xref="RefSeq_mRNA_predicted:XM_533363"
FT                   /db_xref="RefSeq_peptide_predicted:XP_005615333"
FT                   /db_xref="RefSeq_peptide_predicted:XP_022263670"
FT                   /db_xref="RefSeq_peptide_predicted:XP_533363"
FT                   /db_xref="Uniprot/SPTREMBL:E2R204"
FT                   /db_xref="EMBL:AAEX03000011"
FT                   /db_xref="GO:0000398"
FT                   /db_xref="GO:0000398"
FT                   /db_xref="GO:0005634"
FT                   /db_xref="GO:0005682"
FT                   /db_xref="GO:0005829"
FT                   /db_xref="GO:0031965"
FT                   /db_xref="GO:0046540"
FT                   /db_xref="GO:0046540"
FT                   /db_xref="GO:0046540"
FT                   /db_xref="GO:0071005"
FT                   /db_xref="VGNC_trans_name:TXNL4A-201"
FT                   /db_xref="Reactome:R-CFA-72163"
FT                   /db_xref="Reactome:R-CFA-72165"
FT                   /db_xref="Reactome:R-CFA-72172"
FT                   /db_xref="Reactome:R-CFA-72203"
FT                   /db_xref="Reactome:R-CFA-8953854"
FT                   /db_xref="UniParc:UPI0000447A0B"
FT                   /translation="MSYMLPHLHNGWQVDQAILSEEDRVVVIRFGHDWDPTCMKMDEVL
FT                   YSIAEKVKNFAVIYLVDITEVPDFNKMYELYDPCTVMFFFRNKHIMIDLGTGNNNKINW
FT                   AMEDKQEMIDIIETVYRGARKGRGLVVSPKDYSTKYRY"
FT   gene            complement(744461..746178)
FT                   /gene=ENSCAFG00000031133.1
FT                   /locus_tag="HSBP1L1"
FT                   /note="heat shock factor binding protein 1 like 1
FT                   [Source:VGNC Symbol;Acc:VGNC:53725]"
FT   mRNA            join(complement(746112..746178),complement(744461..744552))
FT                   /gene="ENSCAFG00000031133.1"
FT                   /standard_name="ENSCAFT00000045122.1"
FT   CDS             join(complement(746112..746178),complement(744461..744552))
FT                   /gene="ENSCAFG00000031133.1"
FT                   /protein_id="ENSCAFP00000038592.1"
FT                   /note="transcript_id=ENSCAFT00000045122.1"
FT                   /db_xref="RefSeq_mRNA_predicted:XM_003432558"
FT                   /db_xref="RefSeq_peptide_predicted:XP_003432606"
FT                   /db_xref="Uniprot/SPTREMBL:J9NZ72"
FT                   /db_xref="EMBL:AAEX03000011"
FT                   /db_xref="GO:0003714"
FT                   /db_xref="GO:0005634"
FT                   /db_xref="GO:0005737"
FT                   /db_xref="GO:0005829"
FT                   /db_xref="GO:0070370"
FT                   /db_xref="GO:1903507"
FT                   /db_xref="VGNC_trans_name:HSBP1L1-201"
FT                   /db_xref="UniParc:UPI00027479F7"
FT                   /translation="AENLFQELQEHFQALIATLNLRMEEMGSRLEDLQKNVNDLMVQAG
FT                   VEDPVSEQ"
FT   gene            complement(829658..866436)
FT                   /gene=ENSCAFG00000039493.1
FT   misc_RNA        join(complement(865924..866436),complement(846078..846445),
FT                   complement(841686..841787),complement(829658..829667))
FT                   /gene="ENSCAFG00000039493.1"
FT                   /db_xref="RNAcentral:URS0000A92009"
FT                   /note="lincRNA"
FT                   /standard_name="ENSCAFT00000053567.1"
FT   misc_RNA        join(complement(865924..866436),complement(846078..846445),
FT                   complement(843054..843163),complement(841759..841787))
FT                   /gene="ENSCAFG00000039493.1"
FT                   /db_xref="RNAcentral:URS0000A94E9D"
FT                   /note="lincRNA"
FT                   /standard_name="ENSCAFT00000055264.1"
FT   misc_RNA        join(complement(865924..866436),complement(846078..846445),
FT                   complement(845283..845439))
FT                   /gene="ENSCAFG00000039493.1"
FT                   /db_xref="RNAcentral:URS0000AA5288"
FT                   /note="lincRNA"
FT                   /standard_name="ENSCAFT00000058711.1"
FT   gene            complement(886083..886640)
FT                   /gene=ENSCAFG00000028976.1
FT   mRNA            complement(886083..886640)
FT                   /gene="ENSCAFG00000028976.1"
FT                   /standard_name="ENSCAFT00000043602.1"
FT   CDS             complement(886083..886640)
FT                   /gene="ENSCAFG00000028976.1"
FT                   /protein_id="ENSCAFP00000039366.1"
FT                   /note="transcript_id=ENSCAFT00000043602.1"
FT                   /db_xref="Uniprot/SPTREMBL:J9P1E2"
FT                   /db_xref="EMBL:AAEX03000016"
FT                   /db_xref="UniParc:UPI000274763D"
FT                   /translation="MWTGWPMGVPEHCTAPAPYTGRSAQGPSPTSGSAPGPPHTHGPPA
FT                   LGIPPRGPLSTQDYPPTWPPAPRTPLMWAPQQPGPPTQATSTEDHPHATPQHPGLPHPH
FT                   PRGPSAPRTPPCGPSHGSPALGTPPCRPLSTKDPLPPPHPKSYGGWFPGSLFRVLPGPQ
FT                   EDSPPNRAADAQSQHLVAFRCF"


#Example output:

ENSCAFG00000000008	GO:0000398, GO:0000398, GO:0005634, GO:0005682, GO:0005829, GO:0031965, GO:0046540, GO:0046540, GO:0046540, GO:0071005
ENSCAFG00000031133	GO:0003714, GO:0005634, GO:0005737, GO:0005829, GO:0070370, GO:1903507


#command:

$ python Ensembl.dat-to-topGO.db.py -i input.table -o output.tab

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

geneName = "geneName"
dicts = {geneName: []}

outfile = open(args.output, 'w')

with open(args.input) as datafile:
  for line in datafile:
    if line.startswith("FT"):
      words = line.split()
      if '/gene="' in words[1]:
        if dicts[geneName] != []:
            GOprint = ', '.join(str(e) for e in dicts[geneName])
            outfile.write("%s\t%s\n" % (list(dicts.keys())[0], GOprint))
        geneName = words[1].split(".")[0].replace('/gene="', '')
        dicts = {geneName: []}
      elif '/db_xref="GO' in words[1]:
        GO = words[1].replace('/db_xref="', '').replace('"', '')
        dicts[geneName].append(GO)

datafile.close()
outfile.close()

print('Done!')
# dicts = {}
# keys = range(4)
# values = ["Hi", "I", "am", "John"]
# for i in keys:
#         dicts[i] = values[i]
# print(dicts)