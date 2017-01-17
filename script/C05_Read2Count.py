#!/usr/bin/python
import os
import sys
import glob
import string
from Bio import SeqIO
from collections import Counter

def rc(seq):
  seq = str(seq)
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def Processlib(R1file):
  R2file = R1file.replace('_R1_','_R2_')
  R1records = SeqIO.parse(R1file,"fastq")
  R2records = SeqIO.parse(R2file,"fastq")
  R1RandPos = 42
  R2RandPos = 42
  Length    = 18
  muts      = []
  countrecord = 0
  for R1record in R1records:
    countrecord += 1
    if countrecord%100000 == 0:
      print 'Finished Processed %i records' % countrecord
    R2record  = R2records.next()
    R1seq  = R1record.seq
    R2seq  = R2record.seq
    R1roi  = R1seq[R1RandPos:R1RandPos+Length]
    R2roi  = R2seq[R2RandPos:R2RandPos+Length]
    if 'N' in R1roi or 'N' in R2roi: continue
    if R1roi != rc(R2roi): continue
    if R1roi == 'GTTGTCTCCGCTGGTTGG': mut = 'WT'
    elif len(''.join(list(set([R1roi[2],R1roi[5],R1roi[8],R1roi[11],R1roi[14],R1roi[17]]))).replace('G','').replace('T','')) != 0:
      continue
    else: mut =  translation(R1roi)
    muts.append(mut)
  R1records.close()
  R2records.close()
  return Counter(muts)

def Output(MutDict,pop):
  outfile = 'data/'+pop+'.count'
  outfile = open(outfile,'w')
  outfile.write("\t".join(['variant','count'])+"\n")
  for mut in sorted(MutDict.keys(),key=lambda x:MutDict[x],reverse=True):
    outfile.write("\t".join(map(str,[mut,MutDict[mut]]))+"\n")
  outfile.close()

def main():
  R1files = glob.glob('fastq/*_R1_*.fastq')
  MutDict = {}
  sample2pop = {'S1':'Input',
                'S2':'H1R1','S3':'H1R2','S4':'H1R3',
                'S5':'H3R1','S6':'H3R2','S7':'H3R3',
                'S9':'H5R1','S10':'H5R2','S11':'H5R3'}
  for R1file in R1files: 
    pop     = sample2pop[R1file.rsplit('_')[1]]
    MutDict = Processlib(R1file)
    Output(MutDict,pop)
    print 'Finish Processing %s' % pop

if __name__ == "__main__":
  main()
