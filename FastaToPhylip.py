#!/usr/bin/env python

'''
This script converts FASTA to PHYLIP.

The code was inspired by this post https://www.biostars.org/p/56153/#56318

# file.fasta:

>sample1
WYNTNACNCTGSNT---
>sample2
NY-TCACNTCNCNTNNN
>sample3
NNNNNNNNCTNNNNNNN
>sample4
ACCTCACNCYGCNTGGG
>sample5
NCCTCACNCTNCNTGGG
>sample6
NNCTCACNCYNNNTGGG
>sample7
NCCTCACNCTNCNTCCC
>sample8
NNCTCACNCTNNNNGGG 

# file.phy:

 8 17
sample1 WYNTNACNCTGSNT---
sample2 NY-TCACNTCNCNTNNN
sample3 NNNNNNNNCTNNNNNNN
sample4 ACCTCACNCYGCNTGGG
sample5 NCCTCACNCTNCNTGGG
sample6 NNCTCACNCYNNNTGGG
sample7 NCCTCACNCTNCNTCCC
sample8 NNCTCACNCTNNNNGGG

# command:

$ python2 fastaToPhylip.py file.fasta

The output file name will be generated automatically. It will be `file.phy` for this command.

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

# make PHYLIP file name by replacing .fasta with .phy.
phylipName = fastaName.rsplit('.', 1)[0]+'.phy'
if os.path.exists(phylipName):  # check if such PHYLIP exist
  raise Exception('"%s" already exists.' % phylipName)

# create read FASTA generator
sequences = readFasta(fastaName)

# count number of sequences
numSamplesLen = countSampleLength(sequences)

# create read FASTA generator again
sequences = readFasta(fastaName)

# write PHYLIP
phylipFile = open(phylipName, 'w')
phylipFile.write(" %s %s\n" % (numSamplesLen[0], numSamplesLen[1])) # header
for sample in sequences:  # sequences
  phylipFile.write("%s %s\n" % (sample.name, sample.seq))

phylipFile.close()
