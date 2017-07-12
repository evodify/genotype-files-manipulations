#!/usr/bin/env python

'''
This script converts FASTA to tab-delimited table with columns: Chr, Pos, REF.

The code was inspired by this post https://www.biostars.org/p/56153/#56318

# file.fasta:

>sample1
WYNTNACNCTGSNT---
>sample2
NY-TCACNTCNCNTNNN


# file.tab:

Chr Pos REF
sample1 1   W
sample1 2   Y
sample1 3   N
sample1 4   T
sample1 5   N
sample1 6   A
sample1 7   C
sample1 8   N
sample1 9   C
sample1 10  T
sample1 11  G
sample1 12  S
sample1 13  N
sample1 14  T
sample1 15  -
sample1 16  -
sample1 17  -
sample2 1   N
sample2 2   Y
sample2 3   -
sample2 4   T
sample2 5   C
sample2 6   A
sample2 7   C
sample2 8   N
sample2 9   T
sample2 10  C
sample2 11  N
sample2 12  C
sample2 13  N
sample2 14  T
sample2 15  N
sample2 16  N
sample2 17  N

# command:

$ python2 FastaToTab.py file.fasta

The output file name will be generated automatically. It will be `file.tab` for this command.

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

import sys, os, re

class Sequence(object):
  '''Class to store names and sequences'''
  def __init__(self, name, seq):
      self.name = name
      self.seq = seq
  def __len__(self):
      return len(self.seq)

def readFasta(fileName):
  '''Generator to read FASTA file'''
  name = ''
  seq = ''
  with open(fileName) as fastaFile:
    for line in fastaFile:
      line = line.strip()
      if line.startswith('>'):
        if name:
          yield Sequence(name, seq)
        name = line[1:]
        seq = ''
        continue
      seq += line
  yield Sequence(name, seq)

def countSampleLength(SequenceObject):
  for i, sample in enumerate(SequenceObject):
    pass
  return i + 1, len(sample)

# get FASTA file name
fastaName = sys.argv[1]
if not os.path.exists(fastaName):  # check if the file exists
  raise Exception('"%s" does not exists.' % fastaName)

# make Tab file name by replacing .fasta with .tab.
tabName = fastaName.rsplit('.', 1)[0]+'.tab'
if os.path.exists(tabName):  # check if such Tab file exist
  raise Exception('"%s" already exists.' % tabName)

# create read FASTA generator
sequences = readFasta(fastaName)

# count number of sequences
numSamplesLen = countSampleLength(sequences)

# create read FASTA generator again
sequences = readFasta(fastaName)

# write Tab
tabFile = open(tabName, 'w')
tabFile.write("Chr\tPos\tREF\n")
for sample in sequences:
  pos = 0
  for i in sample.seq:
    pos+=1
    tabFile.write("%s\t%s\t%s\n" % (sample.name, pos, i))

tabFile.close()
