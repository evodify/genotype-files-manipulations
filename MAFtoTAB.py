 
#! /usr/bin/env python
'''
Transforms the MAF file to tab file. All gaps will be skipped.

# input:

##maf version=1 scoring=blastz
a score=3680.000000
s scaffold_1     772 292 + 19624517 TCTCTCTCTTCTCTCTTTCATATCCTTTAGAGAGAGAGTCGAACTTTG
s scaffold_5 1072279 229 +  1988289 TCTCTCTCTTCTCTCTATCATATCCTTTAGAGAGAGA---AAGCATTT

# output:

Chr	Pos	REF	TARG
scaffold_1	773	T	T
scaffold_1	774	C	C
scaffold_1	775	T	T
scaffold_1	776	C	C
scaffold_1	777	T	T
scaffold_1	778	C	C
scaffold_1	779	T	T
scaffold_1	780	C	C
scaffold_1	781	T	T
scaffold_1	782	T	T
scaffold_1	783	C	C
scaffold_1	784	T	T
scaffold_1	785	C	C
scaffold_1	786	T	T
scaffold_1	787	C	C
scaffold_1	788	T	T
scaffold_1	789	T	A
scaffold_1	790	T	T
scaffold_1	791	C	C
scaffold_1	792	A	A
scaffold_1	793	T	T
scaffold_1	794	A	A
scaffold_1	795	T	T
scaffold_1	796	C	C
scaffold_1	797	C	C
scaffold_1	798	T	T
scaffold_1	799	T	T
scaffold_1	800	T	T
scaffold_1	801	A	A
scaffold_1	802	G	G
scaffold_1	803	A	A
scaffold_1	804	G	G
scaffold_1	805	A	A
scaffold_1	806	G	G
scaffold_1	807	A	A
scaffold_1	808	G	G
scaffold_1	809	A	A
scaffold_1	813	G	A
scaffold_1	814	A	A
scaffold_1	815	A	G
scaffold_1	816	C	C
scaffold_1	817	T	A
scaffold_1	818	T	T
scaffold_1	819	T	T
scaffold_1	820	G	T

# contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

# command:
python MAFtoTAB.py -i {input} -o {output}

'''

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
args = parser.parse_args()

############################# program #############################

fileoutput = open(args.output, 'w')
fileoutput.write('Chr\tPos\tREF\tTARG\n')

print('Transforming ...')
with open(args.input) as datafile:
  for line in datafile:
    words = line.split()
    if words:
      if (not words[0].startswith("#")):
        param = words[0]
        if param == 's' and linenumber==1:
          REFchr = words[1]
          REFpos = int(words[2])
          REFseq = [i for i in words[6]]
          REFlen = words[3]
          linenumber  = 2
        elif param == 's' and linenumber == 2:
            TARGseq = [i for i in words[6]]
            for i in range(len(REFseq)):
              if REFseq[i] != '-':
                REFpos += 1
                if TARGseq[i] != '-':
                  fileoutput.write('%s\t%s\t%s\t%s\n' % (REFchr, REFpos, REFseq[i], TARGseq[i]))
        elif param == 'a':
          linenumber = 1

datafile.close()
fileoutput.close()
print('Done!')