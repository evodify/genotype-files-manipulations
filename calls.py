#!/usr/bin/python2
'''
This is a python module to operate on call files.

#File examples:

#Two-character code:

CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   1   A   T/A ./. ./. A/A ./. ./. ./. ./.
chr_1   2   C   T/C T/C ./. C/C C/C ./. C/C ./.
chr_1   3   C   C/GCC   C/C ./. C/C C/C C/C C/C C/C
chr_1   4   T   T/T T/T ./. T/T T/T T/T T/T T/T
chr_2   1   A   A/A A/A ./. A/A A/A A/A A/A A/A
chr_2   2   C   C/C C/C ./. C/C C/C C/C C/C C/C
chr_2   3   C   AT/AT   AT/AT   AT/AT   AT/AT   AT/AT   AT/AT   AT/AT   AT/AT
chr_2   4   C   C/C T/T C/C C/C C/C C/C C/C C/C
chr_2   5   T   T/T C/C T/T C/T T/T C/T T/T T/T
chr_3   1   G   G/G ./. ./. G/G ./. ./. ./. ./.
chr_3   2   C   G/C C/C ./. C/C C/C ./. C/C ./.
chr_3   3   CTT CTT/CTT CTT/C   CTT/C   CTT/CTT CTT/CTT CTT/CTT CTT/CTT CTT/CTT
chr_3   4   TA  T/T T/T ./. T/T T/T T/T T/T T/TA
chr_3   5   G   */* G/* ./. G/G G/G G/G C/C G/G


#One-character code:

CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   1   A   W   N   N   A   N   N   N   N
chr_1   2   C   Y   Y   N   C   C   N   C   N
chr_1   3   C   N   C   N   C   C   C   C   C
chr_1   4   T   T   T   N   T   T   T   T   T
chr_2   1   A   A   A   N   A   A   A   A   A
chr_2   2   C   C   C   N   C   C   C   C   C
chr_2   3   C   N   N   N   N   N   N   N   N
chr_2   4   C   C   T   C   C   C   C   C   C
chr_2   5   T   T   C   T   Y   T   Y   T   T
chr_3   1   G   G   N   N   G   N   N   N   N
chr_3   2   C   S   C   N   C   C   N   C   N
chr_3   3   N   N   N   N   N   N   N   N   N
chr_3   4   N   T   T   N   T   T   T   T   N
chr_3   5   G   -   N   N   G   G   G   C   G

Phased:

Note! Some function assumes chromosome number are separated by _.
For example, "chr_1". This format "chr1" may not work all the time.
I am still fixing this issue.

'''

############################# modules #############################

import argparse, sys  # for input options
import collections  # to perform counting
import random  # for randomization

############################# classes  ############################


class CommandLineParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)  # command line syntax errors


class callsParser(object):
    '''
    Parse calls table with genotypes to an object to access
    chromosomes/positions/sequences easily.
    '''
    def __init__(self, filename, samples):
        self.filename = filename
        self.samples = samples
        self.names = []
        self.chrmosomes = []
        self.positions = []
        self.snps = {}

        # read file
        callsFile = open(self.filename, 'r')

        # save samples' names
        header_line = callsFile.readline().split()
        indexS = indexSamples(self.samples, header_line)
        self.names = selectSamples(indexS, header_line)

        # make sequence list
        self.sequences = [[] for i in range(len(self.names))]

        # append sequences
        for line in callsFile:
            words = line.split()
            self.chrmosomes.append(str(words[0]))
            self.positions.append(words[1])
            GT = selectSamples(indexS, words)
            for i in range(len(self.sequences)):
                self.sequences[i].append(GT[i])
            self.snps[words[0]+':'+words[1]] = GT
        callsFile.close()

    def __getitem__(self, i):
        '''
        Enables iteration through chromosomes/positions/sequences by
        index and sample names.
        '''
        if isinstance(i, int):  # if index
            if i > len(self.names) or i < 0:
                raise IndexError("Index is out of range")
            return self.sequences[i]
        else:  # if name
            if i not in self.names:
                raise KeyError("No sequence with name %s", i)
            seqIndex = self.names.index(i)
            return self.sequences[seqIndex]


############################# functions ###########################


def flattenList(complexList):
    '''
    Makes a flat list out of list of lists.
    '''
    flat_list = []
    for sublist in complexList:
        for i in sublist:
            flat_list.append(i)
    return flat_list


def all_missing(genotypes):
    '''
    Check if all genotypes are missing.
    '''
    return all(gt == 'N' for gt in genotypes)


def any_missing(genotypes):
    '''
    Check if any genotype is missing.
    '''
    return any(gt == 'N' for gt in genotypes)


def checkSampleNames(sampleNames, inputFileName):
    '''
    Check if samples names are given and if all sample names are
    present in a header.
    '''
    inputFile = open(inputFileName, 'r')
    inputFile_header = inputFile.readline().split()
    # if no samples specified, use all:
    if sampleNames:
        sampNames = sampleNames.split(',')
        # check if all samples are present in a header
        for sample in sampNames:
            if sample not in inputFile_header:
                raise IOError(
                    'Sample name "%s" is not found in the header' % (sample))
    else:
        sampNames = inputFile_header[2:]
        print 'Sample names are not specified, all will be used ...'
    inputFile.close()
    return sampNames


def indexSamples(sampNames, header_words):
    '''
    extract the index of a given list of sample names.
    '''
    sampIndex = []
    for i in sampNames:
        indnumber = header_words.index(i)
        sampIndex.append(indnumber)
    return sampIndex


def selectSamples(sampIndex, words):
    '''
    extracts column values for given list of indexes.
    '''
    sampWords = []
    for el in sampIndex:
        sampWords.append(words[el])
    return sampWords


def if_all_gt_correct(gt, line):
    '''
    Check if there is any unrecognised genotype.
    '''
    allowed_states = 'AGCTRYMKSWN-*'
    if any(j not in allowed_states for j in gt):
        print('WARNING: unrecognised character in the line -> %s' % line)


def countPerPosition(sampWords, characterToCount):
    '''
    Counts given allele in each position along the genome.
    '''
    count = collections.Counter(sampWords)
    characterCount = count[characterToCount]
    return characterCount


def countHeteroPerPosition(sampWords):
    '''
    Counts heterozygosty in each position along the genome in unphased data.
    '''
    Hcount = 0.0
    for gt in sampWords:
        if gt in "RYSWKM":
            Hcount += 1.0
        elif gt in "ACGTN-":
            continue
        else:
            print('WARNING: character "%s" is not recognized' % gt)
    return Hcount


def if_FixedHetero(sampWords):
    '''
    Returns True if it is a fixed heterozygout (exist in hybrids) and 
    False if not.
    '''
    sampWordsNoNs = []
    for gt in sampWords:  # filter out the missing data
        if gt in "ACGT-RYSWKM":
            sampWordsNoNs.append(gt)
        elif gt == "N":
            continue
        else:
            print('WARNING: character "%s" is not recognized' % gt)
    if all(gt in 'RYMKSW' and gt == sampWordsNoNs[0]
           for gt in sampWordsNoNs):  # check if fixed
        return True
    else:
        return False


def countPerSample(sampWords, countList, characterToCount):
    '''
    Counts Ns (missing data) in each sample.
    '''
    for i in range(len(sampWords)):
        if sampWords[i] == characterToCount:
            countList[i] += 1


def is_polymorphic(sampWords):
    '''
    Check if the set of genotypes is polymorphic.
    '''
    # fist skip missing data
    noNsGT = []
    for i in (sampWords):
        if i != 'N':
            noNsGT.append(i)
    # check if there is polymorphism:
    return any(x in 'RYMKSW' or x != noNsGT[0] for x in noNsGT)


def twoToOne(GT):
    '''
    Converts two character coded genotypes to one character code.
    '''
    GTone = []
    for g in GT:
        if '/' not in g:  # if one character, e.g. the reference column (REF)
            if len(g) != 1:  # if indel
                g = 'N'
        else:  # two character
            if len(g) != 3:  # if indel except single site deletion
                g = 'N'
            elif g[0] == g[2] and g[0] != '.':  #if homozygote and not missing
                if g[0] == '*':  #recode a single site deletion from '*' to '-'
                    g = '-'
                else:
                    g = g[0]
            # single character heterozygouts:
            elif g == 'G/A' or g == 'A/G':
                g = 'R'
            elif g == 'T/C' or g == 'C/T':
                g = 'Y'
            elif g == 'A/C' or g == 'C/A':
                g = 'M'
            elif g == 'G/T' or g == 'T/G':
                g = 'K'
            elif g == 'G/C' or g == 'C/G':
                g = 'S'
            elif g == 'A/T' or g == 'T/A':
                g = 'W'
            else:
                g = 'N'
        GTone.append(g)
    return GTone


def OneToTwo(GT):
    '''
    Converts two character coded genotypes to one character code.
    '''
    GTtwo = []
    for g in GT:
        if '/' not in g:  # if one character, e.g. the reference column (REF)
            if len(g) != 1:  # if indel
                g = './.'
            # single character heterozygouts:
            elif g == 'A':
                g = 'A/A'
            elif g == 'G':
                g = 'G/G'
            elif g == 'C':
                g = 'C/C'
            elif g == 'T':
                g = 'T/T'
            elif g == 'R':
                g = 'A/G'
            elif g == 'Y':
                g = 'C/T'
            elif g == 'M':
                g = 'C/A'
            elif g == 'K':
                g = 'T/G'
            elif g == 'S':
                g = 'C/G'
            elif g == 'W':
                g = 'T/A'
            elif g == 'N':
                g = './.'
            else:
                print(
                    'WARNING: character "%s" is not recognized.'
                    'It will be replaces with N' % g)
                g = './.'
        GTtwo.append(g)
    return GTtwo


def countPositions(fileName):
    '''count number of genomic position in a file'''
    with open(fileName) as f:
        for i in enumerate(f):
            pass
    return i


def processWindow(Chr, FirstPos, LastPos, WindowVal, outputFile):
    '''
    Outputs middle point of a window and a value of this window in a file.
    '''
    #print FirstPos, LastPos
    posP = float(FirstPos) + ((float(LastPos) - float(FirstPos)) / 2.0)
    outputFile.write("%s\t%s\t%s\n" % (Chr, posP, WindowVal))


def pseudoPhase(gt):
    '''
    Randomly splits heterozygouts.
    '''
    # heterozygouts:
    ambR = ['A/G', 'G/A']
    ambY = ['T/C', 'C/T']
    ambM = ['A/C', 'C/A']
    ambK = ['G/T', 'T/G']
    ambS = ['G/C', 'C/G']
    ambW = ['A/T', 'T/A']
    delA = ['A/-', '-/A']
    delT = ['T/-', '-/T']
    delG = ['G/-', '-/G']
    delC = ['C/-', '-/C']
    phasedAlles = []
    for i in gt:
        if i == 'N':
            i = 'N/N'
        elif i == 'A':
            i = 'A/A'
        elif i == 'G':
            i = 'G/G'
        elif i == 'C':
            i = 'C/C'
        elif i == 'T':
            i = 'T/T'
        elif i == '-' or i == '*' or i == "*/*":
            i = '-/-'
        elif i == 'R' or i == 'A/G' or i == 'G/A':
            i = random.choice(ambR)
        elif i == 'Y' or i == 'T/C' or i == 'C/T':
            i = random.choice(ambY)
        elif i == 'M' or i == 'A/C' or i == 'C/A':
            i = random.choice(ambM)
        elif i == 'K' or i == 'G/T' or i == 'T/G':
            i = random.choice(ambK)
        elif i == 'S' or i == 'G/C' or i == 'C/G':
            i = random.choice(ambS)
        elif i == 'W' or i == 'A/T' or i == 'T/A':
            i = random.choice(ambW)
        elif i == 'A/*' or i == '*/A':
            i = random.choice(delA)
        elif i == 'T/*' or i == '*/T':
            i = random.choice(delT)
        elif i == 'G/*' or i == '*/G':
            i = random.choice(delG)
        elif i == 'C/*' or i == '*/C':
            i = random.choice(delC)
        else:
            i = 'N/N'
        phasedAlles.append(i)
    return phasedAlles


def chunks(l, n):
    '''
    Yield successive n-sized chunks from l.
    '''
    for i in xrange(0, len(l), n):
        yield l[i:i + n]


def checkMissing(missingN):
    '''
    Checks if the missing data option is specified.
    Sets no-missing data allowed by default.
    '''
    if missingN:
        allowedNs = missingN
    else:
        allowedNs = 0
        print 'No missing data threshold is specified,'
        'only sites without any missing data will be processed ...'
    return allowedNs


def complementSeq(sequence):
    '''
    Return a complement sequence.
    '''
    complement = []
    for i in sequence:
        if i == 'A' or i == 'a':
            i = 'T'
        elif i == 'T' or i == 't':
            i = 'A'
        elif i == 'G' or i == 'g':
            i = 'C'
        elif i == 'C' or i == 'c':
            i = 'G'
        complement.append(i)
    return complement


def check_if_chr_pos_sorted(callsFile):
    '''
    Check if chr and pos columns are naturally sorted in the calls file.
    '''
    from natsort import natsorted  # to check natural sorting
    print("Checking if %s is naturally sorted..." % callsFile)
    chr = []
    chr_pos = {}
    # load chr and position data
    with open(callsFile, 'r') as f:
        header = f.readline()
        while header.startswith("##"):
            header = f.readline()
        next(f)
        for l in f:
            words = l.rstrip().split()
            chrL = words[0]
            pos = words[1]
            if chrL not in chr:
                chr.append(chrL)
                chr_pos[chrL] = [pos]
            else:
                chr_pos[chrL].append(pos)
    f.close()

    # check if sort order is correct
    if natsorted(chr) != chr:
        raise IOError(
            "Chromosomes are not sorted correctly. "
            "Sort your input with:\nsort -V -k1,1 -k2,2 %s" % callsFile)
    else:
        for k in chr_pos:
            if natsorted(chr_pos[k]) != chr_pos[k]:
                raise IOError(
                    "Genomic positions are not sorted correctly in %s. "\
                    "Sort your input with:\nsort -V -k1,1 -k2,2 %s" % 
                    (k, callsFile))
                break


def convert01toATGC(ref, alt, gt):
    '''
    Convert vcf genotypes (0,1) to nucleotide code (ATGC).
    '''
    altList = alt.split(",")
    reference = [ref] + altList
    gtLetters = reference[int(gt)]
    return gtLetters

def lineCounter(LineNumber):
    LineNumber += 1
    if LineNumber % 100000 == 0:
      print str(LineNumber), "lines processed"
    return LineNumber

def pseudoPhasePED(gt):
  '''
  Randomly splits heterozygouts
  '''
  # heterozygouts:
  ambR = ['A G', 'G A']
  ambY = ['T C', 'C T']
  ambM = ['A C', 'C A']
  ambK = ['G T', 'T G']
  ambS = ['G C', 'C G']
  ambW = ['A T', 'T A']
  delA = ['A -', '- A']
  delT = ['T -', '- T']
  delG = ['G -', '- G']
  delC = ['C -', '- C']
  phasedAlles = []
  for i in gt:
    if i == 'N':
      i = '0 0'
    elif i == 'A':
      i = 'A A'
    elif i == 'G':
      i = 'G G'
    elif i == 'C':
      i = 'C C'
    elif i == 'T':
      i = 'T T'
    elif i == '-' or i == '*' or i == "*/*":
      i = '- -'
    elif i == 'R' or i == 'A/G' or i == 'G/A':
      i = random.choice(ambR)
    elif i == 'Y' or i == 'T/C' or i == 'C/T':
      i = random.choice(ambY)
    elif i == 'M' or i == 'A/C' or i == 'C/A':
      i = random.choice(ambM)
    elif i == 'K' or i == 'G/T' or i == 'T/G':
      i = random.choice(ambK)
    elif i == 'S' or i == 'G/C' or i == 'C/G':
      i = random.choice(ambS)
    elif i == 'W' or i == 'A/T' or i == 'T/A':
      i = random.choice(ambW)
    elif i == 'A/*' or i == '*/A':
      i = random.choice(delA)
    elif i == 'T/*' or i == '*/T':
      i = random.choice(delT)
    elif i == 'G/*' or i == '*/G':
      i = random.choice(delG)
    elif i == 'C/*' or i == '*/C':
      i = random.choice(delC)
    else:
      i = '0 0'
    phasedAlles.append(i)
  return phasedAlles

def phasePED(gt):
  phasedAlles = []
  for i in gt:
    i = i.split("|")
    iP = ' '.join(str(s) for s in i)
    phasedAlles.append(iP)
  return phasedAlles

def familySampleCheck(family_samples, hap_dip_samples):
  ''' To check if samples specified in the options -f, -1n,-2n are the same'''
  if not set(family_samples) == set(hap_dip_samples):
    raise IOError('Sample names in -f differ from those specified in -1n, -2n')

def polarizeGT(genotypes, reference):
    '''
    Polarize a genotype relative to the reference.
    '''
    genotypesPolarized = []
    for gt in genotypes:
        if gt in 'N-.':
            gt = 9
        elif gt == reference:
            gt = 0
        else:
            gt = 1
        genotypesPolarized.append(gt)
    return genotypesPolarized

def is_biallelic(genotypes):
    '''
    Check if a site is biallelic.
    '''
    if all(len(gt) == 1 for gt in genotypes):
        nonMiss = filter(lambda a: a not in 'N.', genotypes)
        for i in range(len(nonMiss)):
            if nonMiss[i] in 'RYMKSW':
                nonMiss[i] = OneToTwo([nonMiss[i]])[0].split('/')
        nonMissGT = flattenList(nonMiss)
        GTset = set(nonMissGT)
        if len(GTset) == 2:
            return True
        else:
            return False
    elif all('/' in gt or '|' in gt for gt in genotypes):
        GTsplit = []
        for gt in genotypes:
            GT = gt.replace('|', '/')
            GT = GT.split('/')
            GTsplit.append(GT)
        GTsplitFlat = flattenList(GTsplit)
        GTsplitFlatNoMiss = filter(lambda a: a not in 'N.', GTsplitFlat)
        GTset = set(GTsplitFlatNoMiss)
        if len(GTset) == 2:
            return True
        else:
            return False
    else:
        raise IOError('Wrong input genotypes %s' % genotypes)
