 
#! /usr/bin/env python
'''
Transforms the MAF file to tab file with Chr Pos of both sequences. All gaps will be skipped.

# input:

##maf version=1 scoring=blastz
a score=3680.000000
s scaffold_1     772 292 + 19624517 TCTCTCTCTTCTCTCTTTCATATCCTTTAGAGAGAGAGTCGAACTTTG
s scaffold_5 1072279 229 +  1988289 TCTCTCTCTTCTCTCTATCATATCCTTTAGAGAGAGA---AAGCATTT

a score=2859.000000
s scaffold_1     166367 101 + 19624517 TTCAAATTTTAAAGATCAAATCTGCAGAAGA----TTTAAAGGTAA--------ATTAAATAAAGACCGTTGCTTT----TTCTTGTTTTCTTCTTCTTCTTCTTAAGTTACAATTT
s MPGU01000319.1 350087 115 -   552602 TTAAAAACGTAAGGTTCAAattca-agaagaagcttttaaAGGTATCatcaattatcaaataaaGACCGttgctttgattttcttactttaacgttctt-ttcaaaatttaaagttt

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
scaffold_1	813	G	A	scaffold_5	1072317
scaffold_1	814	A	A	scaffold_5	1072318
scaffold_1	815	A	G	scaffold_5	1072319
scaffold_1	816	C	C	scaffold_5	1072320
scaffold_1	817	T	A	scaffold_5	1072321
scaffold_1	818	T	T	scaffold_5	1072322
scaffold_1	819	T	T	scaffold_5	1072323
scaffold_1	820	G	T	scaffold_5	1072324
scaffold_1	166368	T	T	MPGU01000319.1	202515
scaffold_1	166369	T	T	MPGU01000319.1	202514
scaffold_1	166370	C	A	MPGU01000319.1	202513
scaffold_1	166371	A	A	MPGU01000319.1	202512
scaffold_1	166372	A	A	MPGU01000319.1	202511
scaffold_1	166373	A	A	MPGU01000319.1	202510
scaffold_1	166374	T	A	MPGU01000319.1	202509
scaffold_1	166375	T	C	MPGU01000319.1	202508
scaffold_1	166376	T	G	MPGU01000319.1	202507
scaffold_1	166377	T	T	MPGU01000319.1	202506
scaffold_1	166378	A	A	MPGU01000319.1	202505
scaffold_1	166379	A	A	MPGU01000319.1	202504
scaffold_1	166380	A	G	MPGU01000319.1	202503
scaffold_1	166381	G	G	MPGU01000319.1	202502
scaffold_1	166382	A	T	MPGU01000319.1	202501
scaffold_1	166383	T	T	MPGU01000319.1	202500
scaffold_1	166384	C	C	MPGU01000319.1	202499
scaffold_1	166385	A	A	MPGU01000319.1	202498
scaffold_1	166386	A	A	MPGU01000319.1	202497
scaffold_1	166387	A	a	MPGU01000319.1	202496
scaffold_1	166388	T	t	MPGU01000319.1	202495
scaffold_1	166389	C	t	MPGU01000319.1	202494
scaffold_1	166390	T	c	MPGU01000319.1	202493
scaffold_1	166391	G	a	MPGU01000319.1	202492
scaffold_1	166393	A	a	MPGU01000319.1	202491
scaffold_1	166394	G	g	MPGU01000319.1	202490
scaffold_1	166395	A	a	MPGU01000319.1	202489
scaffold_1	166396	A	a	MPGU01000319.1	202488
scaffold_1	166397	G	g	MPGU01000319.1	202487
scaffold_1	166398	A	a	MPGU01000319.1	202486
scaffold_1	166398	-	a	MPGU01000319.1	202485
scaffold_1	166398	-	g	MPGU01000319.1	202484
scaffold_1	166398	-	c	MPGU01000319.1	202483
scaffold_1	166398	-	t	MPGU01000319.1	202482
scaffold_1	166399	T	t	MPGU01000319.1	202481
scaffold_1	166400	T	t	MPGU01000319.1	202480
scaffold_1	166401	T	t	MPGU01000319.1	202479
scaffold_1	166402	A	a	MPGU01000319.1	202478
scaffold_1	166403	A	a	MPGU01000319.1	202477
scaffold_1	166404	A	A	MPGU01000319.1	202476
scaffold_1	166405	G	G	MPGU01000319.1	202475
scaffold_1	166406	G	G	MPGU01000319.1	202474
scaffold_1	166407	T	T	MPGU01000319.1	202473
scaffold_1	166408	A	A	MPGU01000319.1	202472
scaffold_1	166409	A	T	MPGU01000319.1	202471
scaffold_1	166409	-	C	MPGU01000319.1	202470
scaffold_1	166409	-	a	MPGU01000319.1	202469
scaffold_1	166409	-	t	MPGU01000319.1	202468
scaffold_1	166409	-	c	MPGU01000319.1	202467
scaffold_1	166409	-	a	MPGU01000319.1	202466
scaffold_1	166409	-	a	MPGU01000319.1	202465
scaffold_1	166409	-	t	MPGU01000319.1	202464
scaffold_1	166409	-	t	MPGU01000319.1	202463
scaffold_1	166410	A	a	MPGU01000319.1	202462
scaffold_1	166411	T	t	MPGU01000319.1	202461
scaffold_1	166412	T	c	MPGU01000319.1	202460
scaffold_1	166413	A	a	MPGU01000319.1	202459
scaffold_1	166414	A	a	MPGU01000319.1	202458
scaffold_1	166415	A	a	MPGU01000319.1	202457
scaffold_1	166416	T	t	MPGU01000319.1	202456
scaffold_1	166417	A	a	MPGU01000319.1	202455
scaffold_1	166418	A	a	MPGU01000319.1	202454
scaffold_1	166419	A	a	MPGU01000319.1	202453
scaffold_1	166420	G	G	MPGU01000319.1	202452
scaffold_1	166421	A	A	MPGU01000319.1	202451
scaffold_1	166422	C	C	MPGU01000319.1	202450
scaffold_1	166423	C	C	MPGU01000319.1	202449
scaffold_1	166424	G	G	MPGU01000319.1	202448
scaffold_1	166425	T	t	MPGU01000319.1	202447
scaffold_1	166426	T	t	MPGU01000319.1	202446
scaffold_1	166427	G	g	MPGU01000319.1	202445
scaffold_1	166428	C	c	MPGU01000319.1	202444
scaffold_1	166429	T	t	MPGU01000319.1	202443
scaffold_1	166430	T	t	MPGU01000319.1	202442
scaffold_1	166431	T	t	MPGU01000319.1	202441
scaffold_1	166431	-	g	MPGU01000319.1	202440
scaffold_1	166431	-	a	MPGU01000319.1	202439
scaffold_1	166431	-	t	MPGU01000319.1	202438
scaffold_1	166431	-	t	MPGU01000319.1	202437
scaffold_1	166432	T	t	MPGU01000319.1	202436
scaffold_1	166433	T	t	MPGU01000319.1	202435
scaffold_1	166434	C	c	MPGU01000319.1	202434
scaffold_1	166435	T	t	MPGU01000319.1	202433
scaffold_1	166436	T	t	MPGU01000319.1	202432
scaffold_1	166437	G	a	MPGU01000319.1	202431
scaffold_1	166438	T	c	MPGU01000319.1	202430
scaffold_1	166439	T	t	MPGU01000319.1	202429
scaffold_1	166440	T	t	MPGU01000319.1	202428
scaffold_1	166441	T	t	MPGU01000319.1	202427
scaffold_1	166442	C	a	MPGU01000319.1	202426
scaffold_1	166443	T	a	MPGU01000319.1	202425
scaffold_1	166444	T	c	MPGU01000319.1	202424
scaffold_1	166445	C	g	MPGU01000319.1	202423
scaffold_1	166446	T	t	MPGU01000319.1	202422
scaffold_1	166447	T	t	MPGU01000319.1	202421
scaffold_1	166448	C	c	MPGU01000319.1	202420
scaffold_1	166449	T	t	MPGU01000319.1	202419
scaffold_1	166450	T	t	MPGU01000319.1	202418
scaffold_1	166452	T	t	MPGU01000319.1	202417
scaffold_1	166453	T	t	MPGU01000319.1	202416
scaffold_1	166454	C	c	MPGU01000319.1	202415
scaffold_1	166455	T	a	MPGU01000319.1	202414
scaffold_1	166456	T	a	MPGU01000319.1	202413
scaffold_1	166457	A	a	MPGU01000319.1	202412
scaffold_1	166458	A	a	MPGU01000319.1	202411
scaffold_1	166459	G	t	MPGU01000319.1	202410
scaffold_1	166460	T	t	MPGU01000319.1	202409
scaffold_1	166461	T	t	MPGU01000319.1	202408
scaffold_1	166462	A	a	MPGU01000319.1	202407
scaffold_1	166463	C	a	MPGU01000319.1	202406
scaffold_1	166464	A	a	MPGU01000319.1	202405
scaffold_1	166465	A	g	MPGU01000319.1	202404
scaffold_1	166466	T	t	MPGU01000319.1	202403
scaffold_1	166467	T	t	MPGU01000319.1	202402
scaffold_1	166468	T	t	MPGU01000319.1	202401

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
          TARGlen = int(words[5])
          TARGseq = [i for i in words[6]]
          TARGstrand = words[4]
          if TARGstrand == "-":
            TARGpos = TARGlen-TARGpos+1
          for i in range(len(REFseq)):
            if REFseq[i] != '-':
              REFpos += 1
            if TARGstrand == "+" and TARGseq[i] != '-':
              TARGpos += 1
            elif TARGstrand == "-" and TARGseq[i] != '-':
              TARGpos -= 1
            if TARGseq[i] != '-':
              fileoutput.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (REFchr, REFpos, REFseq[i], TARGseq[i], TARGchr, TARGpos))
        elif param == 'a':
          linenumber = 1

datafile.close()
fileoutput.close()
print('Done!')