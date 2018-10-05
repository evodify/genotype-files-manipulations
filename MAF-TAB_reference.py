 
#! /usr/bin/env python
'''
Transforms the MAF file to tab file with Chr Pos of both sequences. All gaps will be skipped.

# input:

##maf version=1 scoring=blastz
a score=3680.000000
s scaffold_1     772 292 + 19624517 TCTCTCTCTTCTCTCTTTCATATCCTTTAGAGAGAGAGTCGAACTTTG
s MPGU01000319.1 1072279 229 +  1988289 TCTCTCTCTTCTCTCTATCATATCCTTTAGAGAGAGA---AAGCATTT

a score=2859.000000
s scaffold_1     166367 101 + 19624517 TTCAAATTTTAAAGATCAAATCTGCAGAAGA----TTTAAAGGTAA--------ATTAAATAAAGACCGTTGCTTT----TTCTTGTTTTCTTCTTCTTCTTCTTAAGTTACAATTT
s MPGU01000319.1 350087 115 -   552602 TTAAAAACGTAAGGTTCAAattca-agaagaagcttttaaAGGTATCatcaattatcaaataaaGACCGttgctttgattttcttactttaacgttctt-ttcaaaatttaaagttt

# output:

Chr	Pos	REF	TARG	TARGchr	TARGpos
scaffold_1	5	A	A	MPGU01000658.1	8
scaffold_1	6	A	A	MPGU01000658.1	9
scaffold_1	7	C	C	MPGU01000658.1	10
scaffold_1	8	C	T	MPGU01000658.1	11
scaffold_1	9	C	C	MPGU01000658.1	12
scaffold_1	10	T	C	MPGU01000658.1	13
scaffold_1	11	A	T	MPGU01000658.1	14
scaffold_1	12	A	A	MPGU01000658.1	15
scaffold_1	13	A	A	MPGU01000658.1	16
scaffold_1	14	C	C	MPGU01000658.1	17
scaffold_1	15	C	C	MPGU01000658.1	18
scaffold_1	16	C	C	MPGU01000658.1	19
scaffold_1	17	T	T	MPGU01000658.1	20
scaffold_1	18	A	T	MPGU01000658.1	21
scaffold_1	19	A	A	MPGU01000658.1	22
scaffold_1	20	A	A	MPGU01000658.1	23
scaffold_1	21	C	C	MPGU01000658.1	24
scaffold_1	22	C	C	MPGU01000658.1	25
scaffold_1	23	C	C	MPGU01000658.1	26
scaffold_1	24	T	T	MPGU01000658.1	27
scaffold_1	25	A	A	MPGU01000658.1	28
scaffold_1	26	A	A	MPGU01000658.1	29
scaffold_1	27	A	A	MPGU01000658.1	30
scaffold_1	28	C	C	MPGU01000658.1	31
scaffold_1	29	C	C	MPGU01000658.1	32
scaffold_1	30	C	C	MPGU01000658.1	33
scaffold_1	31	T	C	MPGU01000658.1	34
scaffold_1	32	A	G	MPGU01000658.1	35
scaffold_1	33	A	A	MPGU01000658.1	36
scaffold_1	34	A	A	MPGU01000658.1	37
scaffold_1	38	T	T	MPGU01000658.1	38
scaffold_1	39	A	A	MPGU01000658.1	39
scaffold_1	40	A	A	MPGU01000658.1	40
scaffold_1	41	A	A	MPGU01000658.1	41
scaffold_1	42	C	C	MPGU01000658.1	42
scaffold_1	43	C	A	MPGU01000658.1	43
scaffold_1	44	C	G	MPGU01000658.1	44
scaffold_1	45	T	T	MPGU01000658.1	45
scaffold_1	46	A	A	MPGU01000658.1	46
scaffold_1	47	A	A	MPGU01000658.1	47
scaffold_1	49	C	C	MPGU01000658.1	48
scaffold_1	50	C	C	MPGU01000658.1	49
scaffold_1	51	C	C	MPGU01000658.1	50
scaffold_1	52	T	T	MPGU01000658.1	51
scaffold_1	53	A	A	MPGU01000658.1	52
scaffold_1	54	A	A	MPGU01000658.1	53
scaffold_1	55	A	A	MPGU01000658.1	54
scaffold_1	56	C	C	MPGU01000658.1	55
scaffold_1	57	C	C	MPGU01000658.1	56
scaffold_1	58	C	C	MPGU01000658.1	57
scaffold_1	59	T	T	MPGU01000658.1	58
scaffold_1	60	A	T	MPGU01000658.1	59
scaffold_1	61	A	A	MPGU01000658.1	60
scaffold_1	62	A	A	MPGU01000658.1	61
scaffold_1	63	C	T	MPGU01000658.1	62
scaffold_1	64	C	C	MPGU01000658.1	66
scaffold_1	65	C	C	MPGU01000658.1	67
scaffold_1	66	T	g	MPGU01000658.1	68
scaffold_1	67	A	t	MPGU01000658.1	69
scaffold_1	68	A	a	MPGU01000658.1	70
scaffold_1	69	A	a	MPGU01000658.1	71
scaffold_1	70	C	c	MPGU01000658.1	72
scaffold_1	71	C	c	MPGU01000658.1	73
scaffold_1	72	C	c	MPGU01000658.1	74
scaffold_1	73	T	t	MPGU01000658.1	75
scaffold_1	74	A	a	MPGU01000658.1	76
scaffold_1	75	A	a	MPGU01000658.1	77
scaffold_1	76	A	a	MPGU01000658.1	78
scaffold_1	77	C	c	MPGU01000658.1	79
scaffold_1	78	C	c	MPGU01000658.1	80
scaffold_1	79	C	c	MPGU01000658.1	81
scaffold_1	80	T	t	MPGU01000658.1	82
scaffold_1	81	A	a	MPGU01000658.1	83
scaffold_1	82	A	a	MPGU01000658.1	84
scaffold_1	83	A	a	MPGU01000658.1	85
scaffold_1	84	C	C	MPGU01000658.1	96
scaffold_1	85	C	C	MPGU01000658.1	97
scaffold_1	86	C	C	MPGU01000658.1	98
scaffold_1	87	T	T	MPGU01000658.1	99
scaffold_1	88	A	A	MPGU01000658.1	100
scaffold_1	89	A	A	MPGU01000658.1	101
scaffold_1	90	A	A	MPGU01000658.1	102
scaffold_1	91	C	T	MPGU01000658.1	103
scaffold_1	92	C	A	MPGU01000658.1	104
scaffold_1	93	C	G	MPGU01000658.1	105
scaffold_1	94	T	T	MPGU01000658.1	106
scaffold_1	95	A	A	MPGU01000658.1	107
scaffold_1	96	A	A	MPGU01000658.1	108
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
fileoutput.write('Chr\tPos\tREF\tTARG\tTARGchr\tTARGpos\tstrand\n')

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
              TARGseq[i] = calls.complementSeq(TARGseq[i])[0]
            if TARGseq[i] != '-' and REFseq[i] != '-':
              fileoutput.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (REFchr, REFpos, REFseq[i], TARGseq[i], TARGchr, TARGpos, TARGstrand))
        elif param == 'a':
          linenumber = 1

datafile.close()
fileoutput.close()
print('Done!')