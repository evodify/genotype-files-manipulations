#!/usr/bin/env python
'''
This script extracts lines from a calls file according to scaffold name, start and end positions.
Every interval is output as an separate file.

calls.file:

#CHROM  POS sample1 sample2 sample3 sample4
scaffold_1  113 G   G   G   N
scaffold_1  117 C   C   N   C
scaffold_1  122 C   C   N   C
scaffold_1  137 A   A   T   N
scaffold_1  139 T   T   T   N
scaffold_1  148 A   A   T   N
scaffold_1  161 C   C   C   T
scaffold_1  170 C   T   C   N
scaffold_1  174 A   A   A   N 
scaffold_2  13 G   G   G   T
scaffold_2  17 C   C   C   C
scaffold_2  22 C   C   G   C
scaffold_2  27 A   T   T   N
scaffold_2  29 T   C   T   C
scaffold_2  38 A   A   T   T
scaffold_2  111 C   C   C   T
scaffold_2  140 C   T   C   N
scaffold_2  178 A   A   A   N 

list.file:

#CHROM  start end
scaffold_1 130  150
scaffold_2  10  80

output.file0:
#CHROM  POS sample1 sample2 sample3 sample4
scaffold_1  137 A   A   T   N
scaffold_1  139 T   T   T   N
scaffold_1  148 A   A   T   N

output.file1:
#CHROM  POS sample1 sample2 sample3 sample4
scaffold_2  13 G   G   G   T
scaffold_2  17 C   C   C   C
scaffold_2  22 C   C   G   C
scaffold_2  27 A   T   T   N
scaffold_2  29 T   C   T   C
scaffold_2  38 A   A   T   T


command:

python select_intervals.py -i calls.file -l list.file -o output.file


contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-l', '--list_names', help = 'file containing list of scaffolds and position', type=str, required=True)
parser.add_argument('-t', '--output_type', help = 'whether output should be in one file (one) or separate output for every interval (separate)', type=str, required=False)
args = parser.parse_args()

############################# program #############################
if not args.output_type:
  output_type = "one"
elif args.output_type in ["one", "separate"]: 
  output_type = args.output_type
else:
  raise IOError('The output type should be either "one" or "separate". You specified %s' % args.output_type)

scaflist = open(args.list_names, "r")
header_scaf = scaflist.readline()
IntervalWords = scaflist.readline().split()
IntervalScaf = int(IntervalWords[0].split('_')[1])
IntervalStart = int(IntervalWords[1])
IntervalEnd = int(IntervalWords[2])

if output_type == "separate":
  interval_count = 0
else:
  interval_count = ""

output = open(args.output+str(interval_count), 'w')

counter = 0

with open(args.input) as datafile:
  header = datafile.readline()
  output.write("%s" % header)

  for line in datafile:
    words = line.split()
    scafCalls = int(words[0].split('_')[1])
    posCalls = int(words[1])

    while  scafCalls > IntervalScaf or (scafCalls == IntervalScaf and posCalls > IntervalEnd):
      IntervalWords = scaflist.readline().split()
      if IntervalWords == []:
        break
      else:
        if output_type == "separate":
          interval_count+=1
          output.close()
          output = open(args.output+str(interval_count), 'w')
          output.write("%s" % header)
        IntervalScaf = int(IntervalWords[0].split('_')[1])
        IntervalStart = int(IntervalWords[1])
        IntervalEnd = int(IntervalWords[2])

    if scafCalls == IntervalScaf and posCalls >= IntervalStart and posCalls <= IntervalEnd:
      pass
    else:
      output.write(line)

    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

datafile.close()
output.close()
print('Done!')
