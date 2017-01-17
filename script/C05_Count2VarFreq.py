#!/usr/bin/python
import os
import sys
import glob
import string
from collections import Counter

def hashincounttable(filename):
  MutDict = {}
  infile = open(filename,'r')
  countline = 0
  for line in infile:
    if 'count' in line: continue
    line = line.rstrip().rsplit("\t")
    MutDict[line[0]] = int(line[1])
  infile.close()
  return MutDict

def FreqCal(Dict, mut, Total):
  freq = float(Dict[mut])/float(Total)*1000000 if Dict.has_key(mut) else float(0)
  return freq

def main():
  outfile    = 'data/VariantFreqTable.tsv'
  InputDict  = hashincounttable('data/Input.count')
  H1R1Dict   = hashincounttable('data/H1R1.count')
  H1R2Dict   = hashincounttable('data/H1R2.count')
  H1R3Dict   = hashincounttable('data/H1R3.count')
  H3R1Dict   = hashincounttable('data/H3R1.count')
  H3R2Dict   = hashincounttable('data/H3R2.count')
  H3R3Dict   = hashincounttable('data/H3R3.count')
  H5R1Dict   = hashincounttable('data/H5R1.count')
  H5R2Dict   = hashincounttable('data/H5R2.count')
  H5R3Dict   = hashincounttable('data/H5R3.count')
  muts       = list(set(InputDict.keys()+
                        H1R1Dict.keys()+H1R2Dict.keys()+H1R3Dict.keys()+
                        H3R1Dict.keys()+H3R2Dict.keys()+H3R3Dict.keys()+
                        H5R1Dict.keys()+H5R2Dict.keys()+H5R3Dict.keys()))
  InputTotal = sum(InputDict.values())
  H1R1Total = sum(H1R1Dict.values())
  H1R2Total = sum(H1R2Dict.values())
  H1R3Total = sum(H1R3Dict.values())
  H3R1Total = sum(H3R1Dict.values())
  H3R2Total = sum(H3R2Dict.values())
  H3R3Total = sum(H3R3Dict.values())
  H5R1Total = sum(H5R1Dict.values())
  H5R2Total = sum(H5R2Dict.values())
  H5R3Total = sum(H5R3Dict.values())

  Inputmin   = float(1)/float(InputTotal)*1000000
  outfile       = open(outfile,'w')
  outfile.write("\t".join(['Variant','Input',
                           'H1R1','H1R2','H1R3',
                           'H3R1','H3R2','H3R3',
                           'H5R1','H5R2','H5R3'])+"\n")
  mutcount   = 0
  for mut in muts:
    mutcount += 1
    if mutcount%100000 == 0: print 'Processed %i variants' % mutcount
    Inputfreq = FreqCal(InputDict, mut, InputTotal)
    H1R1freq  = FreqCal(H1R1Dict, mut, H1R1Total)
    H1R2freq  = FreqCal(H1R2Dict, mut, H1R2Total)
    H1R3freq  = FreqCal(H1R3Dict, mut, H1R3Total)
    H3R1freq  = FreqCal(H3R1Dict, mut, H3R1Total)
    H3R2freq  = FreqCal(H3R2Dict, mut, H3R2Total)
    H3R3freq  = FreqCal(H3R3Dict, mut, H3R3Total)
    H5R1freq  = FreqCal(H5R1Dict, mut, H5R1Total)
    H5R2freq  = FreqCal(H5R2Dict, mut, H5R2Total)
    H5R3freq  = FreqCal(H5R3Dict, mut, H5R3Total)
    outfile.write("\t".join(map(str,[mut,Inputfreq,
                                     H1R1freq,H1R2freq,H1R3freq,
                                     H3R1freq,H3R2freq,H3R3freq,
                                     H5R1freq,H5R2freq,H5R3freq]))+"\n")
  outfile.close()

if __name__ == "__main__":
  main()
