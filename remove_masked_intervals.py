#!/usr/bin/env python
'''
This script compares an interval file (BED format) with the interval file of masked regions and remove them.

# test.bed:

Chr Start   End
scaffold_1  1   712
scaffold_1  2790    3155
scaffold_1  10616   11196
scaffold_1  15910   15953
scaffold_1  24547   24981



# masked.bed:

Chr Start   End
scaffold_1  500   700
scaffold_1  10000   20000

# output.bed:

Chr Start   End
scaffold_1  2790    3155
scaffold_1  24547   24981


# command:

python remove_masked_intervals.py -i test.bed -m masked.bed -o output.bed -f 0.0

contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input BED file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-m', '--masked_intervals', help = 'file containing list of coordinates of masked regions in the BED format', type=str, required=True)
parser.add_argument('-f', '--filter', help = 'allowed proportion of masked sites', type=float, required=False)
args = parser.parse_args()

if args.filter == []: # if filter is not specified, skip all overlaps
  filter = 0.0
elif  args.filter >= 0.0 and args.filter <= 1.0:
  filter = args.filter
else: 
  raise IOError('The filter should be of the float type (e.g. 0.1 or 0.7 ). You specified %s' % args.filter)

############################## functions #############################

def getOverlap(a, b):
  return max(0, min(int(a[1]), int(b[1])) - max(int(a[0]), int(b[0])))

############################# program #############################

maskFile = open(args.masked_intervals, "r")
header_scaf = maskFile.readline()
maskWords = maskFile.readline().split()
maskScaf = int(maskWords[0].split('_')[1])
maskCoord = maskWords[1:]
maskStart = int(maskWords[1])
maskEnd = int(maskWords[2])

output = open(args.output, 'w')

counter = 0

with open(args.input) as datafile:
  header = datafile.readline()
  output.write("%s" % header)

  for line in datafile:
    words = line.split()
    inScaff = int(words[0].split('_')[1])
    inCoord = words[1:]
    inStart = int(words[1])
    inEnd = int(words[2])
    overlap = 0

    # read masked coordinates until an overlap
    while  inScaff > maskScaf or (inScaff == maskScaf and inStart > maskEnd):
      maskWords = maskFile.readline().split()
      if maskWords == []:
        break
      else:
        maskScaf = int(maskWords[0].split('_')[1])
        maskCoord = maskWords[1:]
        maskStart = int(maskWords[1])
        maskEnd = int(maskWords[2])

    # find overlap
    if inScaff < maskScaf:
      output.write(line)
      continue
    elif inScaff == maskScaf and maskStart >= inStart and inEnd <= maskEnd:
      overlap += getOverlap(inCoord, maskCoord)
    else:  ## read masked coordinates if they overlap
      while  inScaff == maskScaf and maskStart <= inStart and maskEnd <= inEnd:
        maskWords = maskFile.readline().split()
        if maskWords == []:
          break
        else:
          overlap += getOverlap(inCoord, maskCoord)
          maskScaf = int(maskWords[0].split('_')[1])
          maskCoord = maskWords[1:]
          maskStart = int(maskWords[1])
          maskEnd = int(maskWords[2])

    # check how many sites are masked and omit only good intervals
    overlap += getOverlap(inCoord, maskCoord)
    maskedSize = float(overlap)/float(inEnd - inStart)

    #output.write("%s\t%s\t%s\t%s\n" % (words[0], inStart, inEnd, maskedSize))
    if maskedSize <= filter:
      output.write(line)
      #pass

    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
output.close()
print('Done!')
