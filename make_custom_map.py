#! /usr/bin/python
'''
This script creates a custom map file from a genotype tab/vcf file and a known map.

It is very slow. The loop [for l in range(0,len(maplist)): ...] needs to be replaced with a better approach.

# input.vcf:

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  AM001   AM002   AM127
chr1       29766   .       C       T       .       PASS    .       GT      0|0     0|0     0|0
chr1       29850   .       G       A       .       PASS    .       GT      0|0     0|0     0|0
chr1       30050   .       C       T       .       PASS    .       GT      0|0     0|0     0|0
chr1       30121   .       C       T       .       PASS    .       GT      0|0     0|0     0|0
chr1       30136   .       C       T       .       PASS    .       GT      0|0     0|0     0|0
chr1       30137   .       G       C       .       PASS    .       GT      0|0     0|0     0|0
chr1       30183   .       C       G       .       PASS    .       GT      0|0     0|0     0|0
chr1       30213   .       A       G       .       PASS    .       GT      0|0     0|0     0|0
chr1       30232   .       A       G       .       PASS    .       GT      0|0     0|0     0|0
chr1       30381   .       A       G       .       PASS    .       GT      0|0     0|0     0|0
chr1       30413   .       T       A       .       PASS    .       GT      0|1     0|0     1|0
chr1       30435   .       C       G       .       PASS    .       GT      0|1     0|0     0|0
chr1       30455   .       C       T       .       PASS    .       GT      0|0     0|0     0|0
chr1       30475   .       A       T       .       PASS    .       GT      0|0     0|0     0|0
chr1       30485   .       C       T       .       PASS    .       GT      0|0     0|0     0|0
chr1       30498   .       C       T       .       PASS    .       GT      0|0     0|0     0|0
chr1       30503   .       G       C       .       PASS    .       GT      0|0     0|0     0|0
chr1       30544   .       G       C       .       PASS    .       GT      0|0     0|0     0|0
chr1       30547   .       G       C       .       PASS    .       GT      0|0     0|0     0|0
chr1       30569   .       C       G       .       PASS    .       GT      1|1     1|1     0|1
chr2       215     .       T       G       .       PASS    .       GT      1|0     1|0     1|1
chr2       245     .       A       C       .       PASS    .       GT      0|0     0|0     0|0
chr2       261     .       C       A       .       PASS    .       GT      0|0     0|0     0|0
chr2       280     .       C       A       .       PASS    .       GT      0|0     1|0     0|0
chr2       321     .       G       C       .       PASS    .       GT      0|0     0|0     0|0
chr2       328     .       C       A       .       PASS    .       GT      0|0     0|0     0|0
chr2       345     .       G       A       .       PASS    .       GT      0|0     0|0     0|0
chr2       365     .       A       G       .       PASS    .       GT      0|0     0|0     0|0
chr2       367     .       G       C       .       PASS    .       GT      0|0     0|0     0|0
chr2       373     .       A       C       .       PASS    .       GT      0|0     0|0     0|0
chr2       379     .       C       A       .       PASS    .       GT      0|0     0|0     0|0
chr2       393     .       C       T,G     .       PASS    .       GT      0|0     0|0     0|0
chr2       406     .       C       T       .       PASS    .       GT      0|0     0|0     0|0
chr2       409     .       A       G       .       PASS    .       GT      0|0     0|0     0|0
chr2       418     .       C       T       .       PASS    .       GT      0|0     0|0     0|0
chr2       423     .       G       A       .       PASS    .       GT      0|0     0|0     0|0
chr2       427     .       C       A       .       PASS    .       GT      0|0     0|0     0|0
chr2       445     .       C       T       .       PASS    .       GT      0|0     0|0     0|0
chr2       457     .       G       A       .       PASS    .       GT      0|0     0|0     0|0
chr2       464     .       G       T       .       PASS    .       GT      0|0     0|0     0|0
chr2       467     .       A       C       .       PASS    .       GT      0|0     0|0     0|0
chr2       470     .       G       C       .       PASS    .       GT      0|0     0|0     0|0
chr2       9000    .       G       C       .       PASS    .       GT      0|0     0|0     0|0


# input.map:

POS     POS     RATE(cM/Mb)     MAP(cM)
chr1    30050   3.192200        0.000000
chr1    30157   3.201262        0.000342
chr1    30173   3.212258        0.000393
chr1    30370   3.005506        0.001026
chr1    30406   2.797626        0.001134
chr1    30664   0.785226        0.001856
chr1    32790   0.688866        0.003525
chr1    32974   0.539750        0.003652
chr1    34601   0.008877        0.004530
chr2    296     0.001104        0.000000
chr2    6819    0.001204        0.000007
chr2    8000    0.001404        0.000008
chr3    300   0.192200        0.000000


# output.map:

chr	id	genetic_position	physical_position
chr1	chr1_29766	0.0	29766
chr1	chr1_29850	0.0	29850
chr1	chr1_30050	0.0	30050
chr1	chr1_30121	0.000226934579439	30121
chr1	chr1_30136	0.000274878504673	30136
chr1	chr1_30137	0.000278074766355	30137
chr1	chr1_30183	0.000425131979695	30183
chr1	chr1_30213	0.000521527918782	30213
chr1	chr1_30232	0.000582578680203	30232
chr1	chr1_30381	0.001059	30381
chr1	chr1_30413	0.00115358914729	30413
chr1	chr1_30435	0.00121515503876	30435
chr1	chr1_30455	0.00127112403101	30455
chr1	chr1_30475	0.00132709302326	30475
chr1	chr1_30485	0.00135507751938	30485
chr1	chr1_30498	0.00139145736434	30498
chr1	chr1_30503	0.0014054496124	30503
chr1	chr1_30544	0.00152018604651	30544
chr1	chr1_30547	0.00152858139535	30547
chr1	chr1_30569	0.00159014728682	30569
chr2	chr2_245	0.0	245
chr2	chr2_261	0.0	261
chr2	chr2_280	0.0	280
chr2	chr2_321	2.68281465583e-08	321
chr2	chr2_328	3.43400275947e-08	328
chr2	chr2_345	5.25831672543e-08	345
chr2	chr2_365	7.4045684501e-08	365
chr2	chr2_367	7.61919362257e-08	367
chr2	chr2_373	8.26306913997e-08	373
chr2	chr2_379	8.90694465737e-08	379
chr2	chr2_393	1.04093208646e-07	393
chr2	chr2_406	1.18043844857e-07	406
chr2	chr2_409	1.21263222444e-07	409
chr2	chr2_418	1.30921355205e-07	418
chr2	chr2_423	1.36286984516e-07	423
chr2	chr2_427	1.40579487966e-07	427
chr2	chr2_445	1.59895753488e-07	445
chr2	chr2_457	1.72773263836e-07	457
chr2	chr2_464	1.80285144872e-07	464
chr2	chr2_467	1.83504522459e-07	467
chr2	chr2_470	1.86723900046e-07	470
chr2	chr2_9000	8e-06	9000

# command

$ python make_custom_map.py -i input.vcf -m input.map -o new.map


# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

############################# functions #############################


############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input_vcf', help = 'name of the vcf/tab file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output map file', type=str, required=True)
parser.add_argument('-m', '--input_map', help = 'name of the reference map file', type=str, required=True)
args = parser.parse_args()

############################# program #############################

mapFile = open(args.input_map, 'r')
output = open(args.output, 'w')

mapheader = mapFile.readline()
mapwords = mapFile.readline().split()
mapChrPrev = mapwords[0]
mapPosPrev = int(mapwords[1])
mapMPrev = float(mapwords[3])

mapwords = mapFile.readline().split()
mapChr = mapwords[0]
mapPos = int(mapwords[1])
mapM = float(mapwords[3])

with open(args.input_vcf) as datafile:
  header = datafile.readline()
  output.write("chr\tid\tgenetic_position\tphysical_position\n")

  for line in datafile:
    words = line.split()
    chr = words[0]
    pos = int(words[1])
    # print chr, pos, mapChr, mapPos, mapChrPrev, mapPosPrev
    if chr == mapChrPrev and chr == mapChr:
      if pos < mapPosPrev or pos == mapPosPrev: # before reference SNPs
        #print "before or match", mapPosPrev, pos, mapPos
        output.write("%s\t%s_%s\t%s\t%s\n" % (chr, chr, pos, mapMPrev, pos))
      elif mapPosPrev < pos and pos < mapPos:
        #print "between", mapPosPrev, pos, mapPos
        perBp = (mapM - mapMPrev) / (mapPos - mapPosPrev)
        newM = (pos - mapPosPrev)*perBp + mapMPrev
        output.write("%s\t%s_%s\t%s\t%s\n" % (chr, chr, pos, newM, pos))
      else:
        # print "read", mapPosPrev, pos, mapPos
        try:
          # read new map line
          while chr == mapChr and pos > mapPos:
            mapChrPrev = mapChr
            mapPosPrev = mapPos
            mapMPrev = mapM
            mapwords = mapFile.readline().split()
            mapChr = mapwords[0]
            mapPos = int(mapwords[1])
            mapM = float(mapwords[3])
            #print "have read", mapPosPrev, pos, mapPos
        except:
          pass
        if chr == mapChr and pos == mapPosPrev:
          #print "match", mapPosPrev, pos, mapPos
          output.write("%s\t%s_%s\t%s\t%s\n" % (chr, chr, pos, mapMPrev, pos))
        elif mapPosPrev < pos and pos <= mapPos:
          #print "between", mapPosPrev, pos, mapPos
          perBp = (mapM - mapMPrev) / (mapPos - mapPosPrev)
          newM = (pos - mapPosPrev)*perBp + mapMPrev
          output.write("%s\t%s_%s\t%s\t%s\n" % (chr, chr, pos, newM, pos))
        else: # after reference SNPs
          # print "end", mapChr, mapPos, mapM, mapMPrev
          output.write("%s\t%s_%s\t%s\t%s\n" % (chr, chr, pos, mapMPrev, pos))
          #print "unclassified", mapChrPrev, mapPosPrev,  chr, pos, mapChr, mapPos
    else:
      while chr != mapChrPrev:
        #print "read new line", mapChrPrev, mapPosPrev,  chr, pos, mapChr, mapPos
        mapChrPrev = mapChr
        mapPosPrev = mapPos
        mapMPrev = mapM
        mapwords = mapFile.readline().split()
        mapChr = mapwords[0]
        mapPos = int(mapwords[1])
        mapM = float(mapwords[3])
        #print "have read new line", mapChrPrev, mapPosPrev,  chr, pos, mapChr, mapPos
      if chr == mapChrPrev and chr == mapChr:
        if pos < mapPosPrev or pos == mapPosPrev:  # before reference SNPs
          # print "before or match", mapPosPrev, pos, mapPos
          output.write("%s\t%s_%s\t%s\t%s\n" % (chr, chr, pos, mapMPrev, pos))
        elif mapPosPrev < pos and pos < mapPos:
          # print "between", mapPosPrev, pos, mapPos
          perBp = (mapM - mapMPrev) / (mapPos - mapPosPrev)
          newM = (pos - mapPosPrev) * perBp + mapMPrev
          output.write("%s\t%s_%s\t%s\t%s\n" % (chr, chr, pos, newM, pos))
      else: # end of map chromosome
        output.write("%s\t%s_%s\t%s\t%s\n" % (chr, chr, pos, mapMPrev, pos))

datafile.close()
output.close()
