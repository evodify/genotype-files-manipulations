 
#! /usr/bin/env python
'''
Transforms the MAF file to tab file with Chr Pos of both sequences. All gaps will be skipped.

# input:

##maf version=1 scoring=blastz
a score=3680.000000
s scaffold_1     772 292 + 19624517 TCTCTCTCTTCTCTCTTTCATATCCTTTAGAGAGAGAGTCGAACTTTG
s scaffold_5 1072279 229 +  1988289 TCTCTCTCTTCTCTCTATCATATCCTTTAGAGAGAGA---AAGCATTT

a score=3680.000000
s scaffold_10     7720 24 + 19624517 TTTAGAGAGAGAGTCGAACTTTG
s scaffold_50 10722790 21 -  1988289 TTTAGAGAGAGA---AAGCATTT

# output:

Chr	Pos	REF	TARG	TARGchr	TARGpos
scaffold_1	773	T	T	scaffold_5	1072280
scaffold_1	774	C	C	scaffold_5	1072281
scaffold_1	775	T	T	scaffold_5	1072282
scaffold_1	776	C	C	scaffold_5	1072283
scaffold_1	777	T	T	scaffold_5	1072284
scaffold_1	778	C	C	scaffold_5	1072285
scaffold_1	779	T	T	scaffold_5	1072286
scaffold_1	780	C	C	scaffold_5	1072287
scaffold_1	781	T	T	scaffold_5	1072288
scaffold_1	782	T	T	scaffold_5	1072289
scaffold_1	783	C	C	scaffold_5	1072290
scaffold_1	784	T	T	scaffold_5	1072291
scaffold_1	785	C	C	scaffold_5	1072292
scaffold_1	786	T	T	scaffold_5	1072293
scaffold_1	787	C	C	scaffold_5	1072294
scaffold_1	788	T	T	scaffold_5	1072295
scaffold_1	789	T	A	scaffold_5	1072296
scaffold_1	790	T	T	scaffold_5	1072297
scaffold_1	791	C	C	scaffold_5	1072298
scaffold_1	792	A	A	scaffold_5	1072299
scaffold_1	793	T	T	scaffold_5	1072300
scaffold_1	794	A	A	scaffold_5	1072301
scaffold_1	795	T	T	scaffold_5	1072302
scaffold_1	796	C	C	scaffold_5	1072303
scaffold_1	797	C	C	scaffold_5	1072304
scaffold_1	798	T	T	scaffold_5	1072305
scaffold_1	799	T	T	scaffold_5	1072306
scaffold_1	800	T	T	scaffold_5	1072307
scaffold_1	801	A	A	scaffold_5	1072308
scaffold_1	802	G	G	scaffold_5	1072309
scaffold_1	803	A	A	scaffold_5	1072310
scaffold_1	804	G	G	scaffold_5	1072311
scaffold_1	805	A	A	scaffold_5	1072312
scaffold_1	806	G	G	scaffold_5	1072313
scaffold_1	807	A	A	scaffold_5	1072314
scaffold_1	808	G	G	scaffold_5	1072315
scaffold_1	809	A	A	scaffold_5	1072316
scaffold_1	813	G	A	scaffold_5	1072320
scaffold_1	814	A	A	scaffold_5	1072321
scaffold_1	815	A	G	scaffold_5	1072322
scaffold_1	816	C	C	scaffold_5	1072323
scaffold_1	817	T	A	scaffold_5	1072324
scaffold_1	818	T	T	scaffold_5	1072325
scaffold_1	819	T	T	scaffold_5	1072326
scaffold_1	820	G	T	scaffold_5	1072327
scaffold_10	7721	T	T	scaffold_50	10722789
scaffold_10	7722	T	T	scaffold_50	10722788
scaffold_10	7723	T	T	scaffold_50	10722787
scaffold_10	7724	A	A	scaffold_50	10722786
scaffold_10	7725	G	G	scaffold_50	10722785
scaffold_10	7726	A	A	scaffold_50	10722784
scaffold_10	7727	G	G	scaffold_50	10722783
scaffold_10	7728	A	A	scaffold_50	10722782
scaffold_10	7729	G	G	scaffold_50	10722781
scaffold_10	7730	A	A	scaffold_50	10722780
scaffold_10	7731	G	G	scaffold_50	10722779
scaffold_10	7732	A	A	scaffold_50	10722778
scaffold_10	7736	G	A	scaffold_50	10722774
scaffold_10	7737	A	A	scaffold_50	10722773
scaffold_10	7738	A	G	scaffold_50	10722772
scaffold_10	7739	C	C	scaffold_50	10722771
scaffold_10	7740	T	A	scaffold_50	10722770
scaffold_10	7741	T	T	scaffold_50	10722769
scaffold_10	7742	T	T	scaffold_50	10722768
scaffold_10	7743	G	T	scaffold_50	10722767


# contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

# command:
python MAF-TAB_reference.py -i {input} -o {output}

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
fileoutput.write('Chr\tPos\tREF\tTARG\tTARGchr\tTARGpos\n')

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
          TARGchr = words[1]
          TARGpos = int(words[2])
          TARGseq = [i for i in words[6]]
          TARGstrand = words[4]
          for i in range(len(REFseq)):
            if REFseq[i] != '-':
              REFpos += 1
              if TARGstrand == "+":
                TARGpos += 1
              else:
                TARGpos -= 1
              if TARGseq[i] != '-':
                fileoutput.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (REFchr, REFpos, REFseq[i], TARGseq[i], TARGchr, TARGpos))
        elif param == 'a':
          linenumber = 1

datafile.close()
fileoutput.close()
print('Done!')