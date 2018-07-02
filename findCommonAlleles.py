#!/usr/bin/env python

'''
This script outputs common and rare alleles in a given set of samples.

input file:

CHROM   POS sample1   sample2   sample3   sample4
chr1    1     A/T        A/A       A/T       T/T
chr1    1     C/T        C/C       C/C       C/T     
chr1    1     A/A        A/A       A/T       A/T     
chr1    1     T/T        A/A       A/T       T/T

output file:

#CHR	POS	Common_alleles	Rare_alleles
chr1	1	A,T	N
chr1	1	C	T
chr1	1	A	T
chr1	1	T	A

command:

python findCommonAlleles.py -i test.tab -o out.tab

contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu
'''

############################# modules #############################

import calls # my custom module
import collections

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s', '--sample', help = 'Specify the sample group (optional)', type=str, required=False)
args = parser.parse_args()

# check if samples names are given and if all sample names are present in a header
sNames = calls.checkSampleNames(args.sample, args.input)

############################# program #############################

def all_same(items):
    return all(x == items[0] for x in items[1:])

counter = 0
output = open(args.output, 'w')
output.write("#CHR\tPOS\tCommon_alleles\tRare_alleles\n")

print('Opening the file...')

with open(args.input) as datafile:
    header_words = datafile.readline().split()

    # index samples
    sIndex = calls.indexSamples(sNames, header_words)

    # create lists for output
    sNames = calls.selectSamples(sIndex, header_words)
    
    for line in datafile:
        words = line.split()
        chr_pos = words[0:2]
        chr_posP = '\t'.join(str(e) for e in chr_pos)

        # select samples
        sGT = calls.selectSamples(sIndex, words)

        alleles = sGT
        allelesSplit = [x.split('/') for x in alleles]

        AllallesN = [i for e in allelesSplit for i in e]
        
        if '.' in AllallesN:
            Allalles = [x for x in AllallesN if x != '.']
        else:
            Allalles = AllallesN
            
        numAl = collections.Counter(Allalles)
        numAlV = numAl.values()
        numAlM = numAl.most_common()
        #numAlVS = sorted(numAlV, reverse=True)
        
        if len(Allalles) < 4: # a lot of missing data
            output.write("%s\tN\tN\n" % chr_posP)
        elif all_same(Allalles): # no variation
            output.write("%s\t%s\t-\n" % (chr_posP, Allalles[0]))
        elif all_same(numAlV) : # all alleles with equal frequency
            Alcomm = [i[0] for i in numAlM]
            AlcommP = ','.join(str(e) for e in Alcomm)
            output.write("%s\t%s\tN\n" % (chr_posP, AlcommP))
        elif (len(numAlV)==3) and (numAlM[0][1]==numAlM[1][1]): # two alleles with equal frequency
            Alcomm = [i[0] for i in numAlM[0:2]]
            AlcommP = ','.join(str(e) for e in Alcomm)
            output.write("%s\t%s\t%s\n" % (chr_posP, AlcommP, numAlM[2][0]))
        elif (len(numAlV)==4) and (numAlM[0][1]==numAlM[1][1]) and (numAlM[0][1]!=numAlM[2][1]): # three alleles with equal frequency
            Alcomm = [i[0] for i in numAlM[0:2]]
            AlcommP = ','.join(str(e) for e in Alcomm)
            Alrare = [i[0] for i in numAlM[2:]]
            AlrareP = ','.join(str(e) for e in Alrare)
            output.write("%s\t%s\t%s\n" % (chr_posP, AlcommP, AlrareP))
        elif (len(numAlV)==4) and (numAlM[0][1]==numAlM[1][1]) and (numAlM[0][1]==numAlM[2][1]): # three alleles with equal frequency
            Alcomm = [i[0] for i in numAlM[0:3]]
            AlcommP = ','.join(str(e) for e in Alcomm)
            output.write("%s\t%s\t%s\n" % (chr_posP, AlcommP, numAlM[3][0]))
        else: # one common allele
            Alrare = [i[0] for i in numAlM[1:]]
            AlrareP = ','.join(str(e) for e in Alrare)
            output.write("%s\t%s\t%s\n" % (chr_posP, numAlM[0][0], AlrareP))
    
        counter += 1
        if counter % 1000000 == 0:
            print str(counter), "lines processed"
        
datafile.close()
output.close()
print('Done!')
