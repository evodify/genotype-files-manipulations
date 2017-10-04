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


Note! Chromosome number must be separated by _.
For example, chr_1 - correct, chr1 - incorrect.

'''

############################# modules #############################

import argparse, sys # for input options
import collections # to perform counting
import random # for randomization

############################# classes  ############################

class CommandLineParser(argparse.ArgumentParser): 
   def error(self, message):
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

############################# functions ###########################

def all_missing(genotypes):
  '''check if all genotypes are missing'''
  return all(gt == 'N' for gt in genotypes)


def any_missing(genotypes):
  '''check if any genotype is missing'''
  return any(gt == 'N' for gt in genotypes)


def checkSampleNames(sampleNames, inputFileName):
  '''check if samples names are given and if all sample names are present in a header'''
  inputFile = open(inputFileName, 'r')
  inputFile_header = inputFile.readline().split()
  # if no samples specified, use all:
  if sampleNames:  
    sampNames = sampleNames.split(',')
    # check if all samples are present in a header
    for sample in sampNames:
      if sample not in inputFile_header:
        raise IOError('Sample name "%s" is not found in the header' % (sample))
  else:
    sampNames = inputFile_header[2:]
    print 'Sample names are not specified, all will be used ...'
  inputFile.close()
  return sampNames


def indexSamples(sampNames, header_words):
  ''' extract the index of a given list of sample names'''
  sampIndex = []
  for i in sampNames:
    indnumber = header_words.index(i)
    sampIndex.append(indnumber)
  return sampIndex


def selectSamples(sampIndex, words):
  '''extracts column values for given list of indexes'''
  sampWords = []
  for el in sampIndex:
    sampWords.append(words[el])
  return sampWords

def if_all_gt_correct(gt, line):
  '''Check if there is any unrecognised genotype '''
  allowed_states = 'AGCTRYMKSWN-*'
  if any(j not in allowed_states for j in gt):
    print('WARNING: unrecognised character in the line -> %s' % line)


def countPerPosition(sampWords, characterToCount):
  '''Counts given allele in each position along the genome '''
  count = collections.Counter(sampWords)
  characterCount = count[characterToCount]
  return characterCount


def countHeteroPerPosition(sampWords):
  '''Counts heterozygosty in each position along the genome in unphased data '''
  Hcount = 0.0
  for gt in sampWords:
    if gt in "RYSWKM":
      Hcount += 1.0
    elif gt in "ACGTN-":
      continue
    else:
       print('WARNING: character "%s" is not recognized' % gt)
  return Hcount


def countPerSample(sampWords, countList, characterToCount):
  '''Counts Ns (missing data) in each sample'''
  for i in range(len(sampWords)):
    if sampWords[i] == characterToCount:
      countList[i] += 1


def is_polymorphic(sampWords):
  ''' check if the set of genotypes is polymorphic '''
  # fist skip missing data
  noNsGT = []
  for i in (sampWords):
    if i != 'N':
      noNsGT.append(i)
  # check if there is polymorphism:
  return any(x in 'RYMKSW' or x != noNsGT[0] for x in noNsGT)


def twoToOne(GT):
  ''' converts two character coded genotypes to one character code '''
  GTone = []
  for g in GT:
    if '/' not in g:  # if one character, e.g. the reference column (REF)
      if len(g) != 1:  # if indel
        g = 'N'
    else:  # two character
      if len(g) != 3:  # if indel except single site deletion 
        g = 'N'
      elif g[0] == g[2] and g[0] != '.': # if homozygote and not missing data (./.)
        if g[0] == '*':  # recode a single site deletion from '*' to '-'
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
        g =  'K'
      elif g == 'G/C' or g == 'C/G':
        g = 'S'
      elif g == 'A/T' or g == 'T/A':
        g = 'W'
      else:
        g = 'N'
    GTone.append(g)
  return GTone


def countPositions(fileName):
  '''count number of genomic position in a file'''
  with open(fileName) as f:
      for i, l in enumerate(f):
          pass
  return i


def processWindow(Chr, FirstPos, LastPos, WindowVal, outputFile):
  ''' Outputs middle point of a window and a value of this window in a file  '''
  #print FirstPos, LastPos
  posP = float(FirstPos)+((float(LastPos)-float(FirstPos))/2.0)
  outputFile.write("%s\t%s\t%s\n" % (Chr, posP, WindowVal))


def pseufoPhase(gt):
  ''' Randomly splits heterozygouts '''
  # heterozygouts:
  ambR = ['A\tG', 'G\tA']
  ambY = ['T\tC', 'C\tT']
  ambM = ['A\tC', 'C\tA']
  ambK = ['G\tT', 'T\tG']
  ambS = ['G\tC', 'C\tG']
  ambW = ['A\tT', 'T\tA']
  phasedAlles = []
  for i in gt:
    if i == 'N':
      i = 'N\tN'
    elif i == 'A':
      i = 'A\tA'
    elif i == 'G':
      i = 'G\tG'
    elif i == 'C':
      i = 'C\tC'
    elif i == 'T':
      i = 'T\tT'
    elif i == '-' or i == '*':
      i = '-\t-'
    elif i == 'R':
      i = random.choice(ambR)
    elif i == 'Y':
      i = random.choice(ambY)
    elif i == 'M':
      i = random.choice(ambM)
    elif i == 'K':
      i = random.choice(ambK)
    elif i == 'S':
      i = random.choice(ambS)
    elif i == 'W':
      i = random.choice(ambW)
    phasedAlles.append(i)
  return phasedAlles