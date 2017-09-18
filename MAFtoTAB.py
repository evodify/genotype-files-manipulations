 
#! /usr/bin/env python
'''
Transforms the MAF file to tab file. All gaps will be skipped

# input:

##maf version=1 scoring=blastz
a score=3680.000000
s CR_scaffold_1     772 292 + 19624517 TCTCTCTCTTCTCTCTTTCATATCCTTTAGAGAGAGAGTCGAACTTTGTTCTTTAGTCCGCCGTCCGGAAACTCCGGCAAGGCCGTCCGGTTAGTTTTAACCGGC-TTCAGCCCAACATCCATTGGCAATTCTCATGCTATGGTTGGTTGGTTCTCTTCTTGACATCTTCTGACTCTGTTAGCGGAGATGGTCGG-TTCTCATCCGGCAAAT-CGACTTCTACATGCGGT-TTGGCCGGCGTAATAGTTGCTTTTGTTTTGAATCGGTCTTAAACTGATTGTTATAGTTTGTTTTT
s NP_scaffold_5 1072279 229 +  1988289 TCTCTCTCTTCTCTCTATCATATCCTTTAGAGAGAGA---AAGCATTTTCTTTTGTTTCGCCG---------------------ATCAGTTTATTTTCGGTGAATATGCATCCCAATA--------------------------TAGTTGGGCTTTTGCCGAGCATCTTCTG---------TTGAAGATGGCTGGCTTCTCTTCCGGCAACTCCGGCTTATTTCAGCGGCGTTGGCCGGCTTAGCGTCTG--------GTAACTCCGGCTTACATCAGCGGCGTTAGCTGGCTTCT

# output:

Chr     Pos     REF     TARG
scaffold        773     T       T
scaffold        774     C       C
scaffold        775     T       T
scaffold        776     C       C
scaffold        777     T       T
scaffold        778     C       C
scaffold        779     T       T
scaffold        780     C       C
scaffold        781     T       T
scaffold        782     T       T
scaffold        783     C       C
scaffold        784     T       T
scaffold        785     C       C
scaffold        786     T       T
scaffold        787     C       C
scaffold        788     T       T
scaffold        789     T       A
scaffold        790     T       T
scaffold        791     C       C
scaffold        792     A       A
scaffold        793     T       T
scaffold        794     A       A
scaffold        795     T       T
scaffold        796     C       C
scaffold        797     C       C
scaffold        798     T       T
scaffold        799     T       T
scaffold        800     T       T
scaffold        801     A       A
scaffold        802     G       G
scaffold        803     A       A
scaffold        804     G       G
scaffold        805     A       A
scaffold        806     G       G
scaffold        807     A       A
scaffold        808     G       G
scaffold        809     A       A
scaffold        813     G       A
scaffold        814     A       A
scaffold        815     A       G
scaffold        816     C       C
scaffold        817     T       A
scaffold        818     T       T
scaffold        819     T       T
scaffold        820     G       T
scaffold        821     T       T
scaffold        822     T       C
scaffold        823     C       T
scaffold        824     T       T
scaffold        825     T       T
scaffold        826     T       T
scaffold        827     A       G
scaffold        828     G       T
scaffold        829     T       T
scaffold        830     C       T
scaffold        831     C       C
scaffold        832     G       G
scaffold        833     C       C
scaffold        834     C       C
scaffold        835     G       G
scaffold        857     G       A
scaffold        858     T       T
scaffold        859     C       C
scaffold        860     C       A
scaffold        861     G       G
scaffold        862     G       T
scaffold        863     T       T
scaffold        864     T       T
scaffold        865     A       A
scaffold        866     G       T
scaffold        867     T       T
scaffold        868     T       T
scaffold        869     T       T
scaffold        870     T       C
scaffold        871     A       G
scaffold        872     A       G
scaffold        873     C       T
scaffold        874     C       G
scaffold        875     G       A
scaffold        876     G       A
scaffold        877     C       T
scaffold        878     T       T
scaffold        879     T       G
scaffold        880     C       C
scaffold        881     A       A
scaffold        882     G       T
scaffold        883     C       C
scaffold        884     C       C
scaffold        885     C       C
scaffold        886     A       A
scaffold        887     A       A
scaffold        888     C       T
scaffold        889     A       A
scaffold        916     T       T
scaffold        917     G       A
scaffold        918     G       G
scaffold        919     T       T
scaffold        920     T       T
scaffold        921     G       G
scaffold        922     G       G
scaffold        923     T       G
scaffold        924     T       C
scaffold        925     C       T
scaffold        926     T       T
scaffold        927     C       T
scaffold        928     T       T
scaffold        929     T       G
scaffold        930     C       C
scaffold        931     T       C
scaffold        932     T       G
scaffold        933     G       A
scaffold        934     A       G
scaffold        935     C       C
scaffold        936     A       A
scaffold        937     T       T
scaffold        938     C       C
scaffold        939     T       T
scaffold        940     T       T
scaffold        941     C       C
scaffold        942     T       T
scaffold        943     G       G
scaffold        953     G       T
scaffold        954     C       T
scaffold        955     G       G
scaffold        956     G       A
scaffold        957     A       A
scaffold        958     G       G
scaffold        959     A       A
scaffold        960     T       T
scaffold        961     G       G
scaffold        962     G       G
scaffold        963     T       C
scaffold        964     C       T
scaffold        965     G       G
scaffold        966     G       G
scaffold        967     T       T
scaffold        968     T       T
scaffold        969     C       C
scaffold        970     T       T
scaffold        971     C       C
scaffold        972     A       T
scaffold        973     T       T
scaffold        974     C       C
scaffold        975     C       C
scaffold        976     G       G
scaffold        977     G       G
scaffold        978     C       C
scaffold        979     A       A
scaffold        980     A       A
scaffold        981     A       C
scaffold        982     T       T
scaffold        983     C       C
scaffold        984     G       G
scaffold        985     A       G
scaffold        986     C       C
scaffold        987     T       T
scaffold        988     T       T
scaffold        989     C       A
scaffold        990     T       T
scaffold        991     A       T
scaffold        992     C       T
scaffold        993     A       C
scaffold        994     T       A
scaffold        995     G       G
scaffold        996     C       C
scaffold        997     G       G
scaffold        998     G       G
scaffold        999     T       C
scaffold        1000    T       T
scaffold        1001    T       T
scaffold        1002    G       G
scaffold        1003    G       G
scaffold        1004    C       C
scaffold        1005    C       C
scaffold        1006    G       G
scaffold        1007    G       G
scaffold        1008    C       C
scaffold        1009    G       T
scaffold        1010    T       T
scaffold        1011    A       A
scaffold        1012    A       G
scaffold        1013    T       C
scaffold        1014    A       G
scaffold        1015    G       T
scaffold        1016    T       C
scaffold        1017    T       T
scaffold        1018    G       G
scaffold        1027    T       G
scaffold        1028    T       T
scaffold        1029    G       A
scaffold        1030    A       A
scaffold        1031    A       C
scaffold        1032    T       T
scaffold        1033    C       C
scaffold        1034    G       C
scaffold        1035    G       G
scaffold        1036    T       G
scaffold        1037    C       C
scaffold        1038    T       T
scaffold        1039    T       T
scaffold        1040    A       A
scaffold        1041    A       C
scaffold        1042    A       A
scaffold        1043    C       T
scaffold        1044    T       C
scaffold        1045    G       A
scaffold        1046    A       G
scaffold        1047    T       C
scaffold        1048    T       G
scaffold        1049    G       G
scaffold        1050    T       C
scaffold        1051    T       G
scaffold        1052    A       T
scaffold        1053    T       T
scaffold        1054    A       A
scaffold        1055    G       G
scaffold        1056    T       C
scaffold        1057    T       T
scaffold        1058    T       G
scaffold        1059    G       G
scaffold        1060    T       C
scaffold        1061    T       T
scaffold        1062    T       T
scaffold        1063    T       C
scaffold        1064    T       T

# contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

# command:
python MAFtoTAB.py -i {input} -o {output}

'''

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
args = parser.parse_args()

############################# program #############################

fileoutput = open(args.output, 'w')
fileoutput.write('Chr\tPos\tREF\tTARG\n')

print('Transforming ...')
with open(args.input) as datafile:
  for line in datafile:
    if not line.strip():
      continue
    else:
      words = line.split()
      param = words[0]
    if param == 's':
      CH = words[1].split('_')
      spCH = CH[0]
      asCH = CH[1]
      if spCH == 'CR':
        CRseq = [i for i in words[6]]
        CRpos = int(words[2])
        CRlen = words[3]
        CHR = asCH
      else:
        NPseq = [i for i in words[6]]
        NPpos = words[2]
        NPlen = words[3]
        for i in range(len(CRseq)):
          if CRseq[i] != '-':
            CRpos += 1
            if NPseq[i] != '-':
              fileoutput.write('%s\t%s\t%s\t%s\n' % (CHR, CRpos, CRseq[i], NPseq[i]))            

datafile.close()
fileoutput.close()
print('Done!')