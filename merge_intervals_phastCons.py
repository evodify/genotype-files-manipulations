#!/usr/bin/env python
'''
Annotates given interval list with phastCons scrores.

interval.file:

CHROM   POS     END     Wolf_Arc_GLM
chr6    1       1000    0.673519481672333
chr6    1001    2000    0.202736151000017
chr6    2001    3000    0.176488191035573
chr6    3001    4000    0.461673796783709
chr6    4001    5000    0.846457671308096
chr6    5001    6000    0.000172147557810159
chr6    6001    7000    0.459624545847323
chr6    7001    8000    0.461650422787001
chr6    8001    9000    0.946297570672762

phastCons.bed:

chr6   311 324  lod=16  268
chr6   2061 2066  lod=13  247
chr6   2100 2206  lod=376 580
chr6   2653 2728  lod=19  285
chr6   2908 3028  lod=348 573
chr6   3549 3788  lod=900 667
chr6   3589 3784  lod=747 648
chr6   3790 3797  lod=18  280
chr6   5986 6012  lod=34  343
chr6   6019 6039  lod=31  333

output.file:

CHROM	POS	END	Wolf_Arc_GLM	mean_lod	max_lod	min_lod	mean_cons	max_cons	min_cons
chr6	1	1000	0.673519481672333	16.0	16.0	16.0	268.0	268.0	268.0
chr6	1001	2000	0.202736151000017	NA	NA	NA	NA	NA	NA
chr6	2001	3000	0.176488191035573	189.0	376.0	13.0	421.25	580.0	247.0
chr6	3001	4000	0.461673796783709	503.25	900.0	18.0	542.0	667.0	280.0
chr6	4001	5000	0.846457671308096	NA	NA	NA	NA	NA	NA
chr6	5001	6000	0.000172147557810159	34.0	34.0	34.0	343.0	343.0	343.0
chr6	6001	7000	0.459624545847323	32.5	34.0	31.0	338.0	343.0	333.0
chr6	7001	8000	0.461650422787001	NA	NA	NA	NA	NA	NA
chr6	8001	9000	0.946297570672762	NA	NA	NA	NA	NA	NA

command:

$ python merge_intervals-phastCons.py \
    -r interval.file \
    -a phastCons.bed \
    -o output.file

contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls  # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument(
    '-i',
    '--input',
    help='name of the reference input file',
    type=str,
    required=True)
parser.add_argument(
    '-a',
    '--annotation',
    help='name of the annotation file to introduce into the reference input',
    type=str,
    required=True)
parser.add_argument(
    '-o',
    '--output',
    help='name of the output file',
    type=str, required=True)
args = parser.parse_args()

############################# program #############################


with open(args.annotation) as annotsFile:
  annotsWords = annotsFile.readline().rstrip().split("\t")
  InterDicts = {annotsWords[0]: [annotsWords[1:]]}
  prevChr = annotsWords[0]
  for gline in annotsFile:
    glineWords = gline.rstrip().split("\t")
    chr = glineWords[0]
    if chr != prevChr:
      InterDicts[chr] = [glineWords[1:]]
    else:
      InterDicts[chr].append(glineWords[1:])
    prevChr = chr
annotsFile.close()

lodList = []
consList = []
prevIntervalScaf = 'NA'

with open(args.input) as intervalFile:
  header_words = intervalFile.readline()
  output = open(args.output, 'w')
  output.write("%s\tmean_lod\tmax_lod_phastCons\tmin_lod_phastCons\t"
               "mean_cons_phastCons\tmax_cons_phastCons\tmin_cons_phastCons\n" 
                % header_words.rstrip())
  for line in intervalFile:
    intervalWords = line.split()
    intervalScaf = intervalWords[0]
    intervalStart = int(intervalWords[1])
    intervalEnd = int(intervalWords[2])
    if intervalScaf != prevIntervalScaf:
      try:
        chrInterDict = InterDicts[intervalScaf]
      except:
        intervalWordsP = '\t'.join(str(e) for e in intervalWords)
        output.write('%s\tNA\n' % (intervalWordsP))
        continue
    prevIntervalScaf = intervalScaf

    for annotCoord in chrInterDict:
      annotStart = int(annotCoord[0])
      annotEnd = int(annotCoord[1])
      if (intervalStart <= annotStart and annotStart <= intervalEnd) or \
         (intervalStart <= annotEnd and  annotEnd <= intervalEnd) or \
         (annotStart <= intervalStart and intervalEnd <= annotEnd):
        lodList.append(float(annotCoord[2].replace('lod=', '')))
        consList.append(float(annotCoord[3]))

    if intervalWords != []:
      intervalWordsP = '\t'.join(str(e) for e in intervalWords)
      if lodList==[]:
        output.write('%s\tNA\tNA\tNA\tNA\tNA\tNA\n' % intervalWordsP)
      else:
        lodListMean = sum(lodList)/len(lodList)
        lodListMax = max(lodList)
        lodListMin = min(lodList)
        consListMean = sum(consList)/len(consList)
        consListMax = max(consList)
        consListMin = min(consList)
        output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
                    (intervalWordsP, lodListMean, lodListMax, lodListMin,
                     consListMean, consListMax, consListMin))
    lodList = []
    consList = []

intervalFile.close()
output.close()
print('Done!')
