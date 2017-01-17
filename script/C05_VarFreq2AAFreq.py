#!/usr/bin/python
import os
import sys
import glob

def aafrequency(filename,mutaas,pops):
  infile = open(filename,'r')
  FreqDict = {}
  for pop in pops: 
    FreqDict[pop] = {}
    for mutaa in mutaas:
      FreqDict[pop][mutaa] = 0
  for line in infile.xreadlines():
    if 'Variant' in line: continue
    line = line.rstrip().rsplit("\t")
    mut  = line[0]
    if mut == 'WT': continue
    freq_input = float(line[1])
    freq_H1R1  = float(line[2])
    freq_H1R2  = float(line[3])
    freq_H1R3  = float(line[4])
    freq_H3R1  = float(line[5])
    freq_H3R2  = float(line[6])
    freq_H3R3  = float(line[7])
    freq_H5R1  = float(line[8])
    freq_H5R2  = float(line[9])
    freq_H5R3  = float(line[10])
    for p in range(len(mut)):
      aa   = mut[p]
      resi = str(p+1)
      ID   = resi+aa
      if ID not in mutaas: continue
      FreqDict['Input'][ID] += freq_input
      FreqDict['H1R1'][ID] += freq_H1R1
      FreqDict['H1R2'][ID] += freq_H1R2
      FreqDict['H1R3'][ID] += freq_H1R3
      FreqDict['H3R1'][ID] += freq_H3R1
      FreqDict['H3R2'][ID] += freq_H3R2
      FreqDict['H3R3'][ID] += freq_H3R3
      FreqDict['H5R1'][ID] += freq_H5R1
      FreqDict['H5R2'][ID] += freq_H5R2
      FreqDict['H5R3'][ID] += freq_H5R3
  infile.close()
  return FreqDict

def genmutaa():
  aas       = ['E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','_']
  residues  = [1,2,3,4,5,6]
  mutaas    = []
  for residue in residues:
    for aa in aas:
      mutaas.append(str(residue)+aa)
  return mutaas

def Renormalize(FreqDict):
  for pop in FreqDict.keys():
    rescalefactor = float(sum(FreqDict[pop].values()))/6
    for ID in FreqDict[pop].keys():
      FreqDict[pop][ID] = FreqDict[pop][ID]/rescalefactor
  return FreqDict

def main():
  filename = 'data/VariantFreqTable.tsv'
  outfile  = 'data/AAFreqTable.tsv'
  mutaas   = genmutaa()
  mutaas   = list(set(mutaas)-set(['6E','6D','6R','6K','6H','6Q','6N','6S','6T','6P','6G','6A','6V','6I','6M']))
  pops     = ['Input','H1R1','H1R2','H1R3','H3R1','H3R2','H3R3','H5R1','H5R2','H5R3']
  FreqDict = aafrequency(filename,mutaas,pops)
  FreqDict = Renormalize(FreqDict)
  outfile  = open(outfile,'w')
  outfile.write('pos'+"\t"+'aa'+"\t"+"\t".join(pops)+"\n")
  freqmaxvalue   = 0
  for mutaa in mutaas:
    out = [mutaa[0],mutaa[1]]
    for pop in pops:
      freq   = FreqDict[pop][mutaa]
      out.append(str(freq))
      if freq   > freqmaxvalue:   freqmaxvalue = freq
    outfile.write("\t".join(out)+"\n")
  outfile.close()
  print "Largest Frequency Value = %f" % freqmaxvalue

if __name__ == "__main__":
  main()
