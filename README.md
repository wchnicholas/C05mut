## YEAST DISPLAY SCREENING FOR C05 VARIANTS
### FILES
* All sequencing raw reads, which can be downloaded from NIH SRA database [PRJNA326694](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA326694), should placed in fastq/ folder:
  * Input library: Wu-1\_S1\_L001\_R1\_001.fastq and Wu-1\_S1\_L001\_R2\_001.fastq
  * H1 Round 1 library: Wu-1\_S2\_L001\_R1\_001.fastq and Wu-1\_S2\_L001\_R2\_001.fastq
  * H1 Round 2 library: Wu-1\_S3\_L001\_R1\_001.fastq and Wu-1\_S3\_L001\_R2\_001.fastq
  * H1 Round 3 library: Wu-1\_S4\_L001\_R1\_001.fastq and Wu-1\_S4\_L001\_R2\_001.fastq
  * H3 Round 1 library: Wu-1\_S5\_L001\_R1\_001.fastq and Wu-1\_S5\_L001\_R2\_001.fastq
  * H3 Round 2 library: Wu-1\_S6\_L001\_R1\_001.fastq and Wu-1\_S6\_L001\_R2\_001.fastq
  * H3 Round 3 library: Wu-1\_S7\_L001\_R1\_001.fastq and Wu-1\_S7\_L001\_R2\_001.fastq
  * H5 Round 1 library: Wu-1\_S9\_L001\_R1\_001.fastq and Wu-1\_S9\_L001\_R2\_001.fastq
  * H5 Round 2 library: Wu-1\_S10\_L001\_R1\_001.fastq and Wu-1\_S10\_L001\_R2\_001.fastq
  * H5 Round 3 library: Wu-1\_S11\_L001\_R1\_001.fastq and Wu-1\_S11\_L001\_R2\_001.fastq
* data/\*count: Counting for each variant in each library
* data/VariantFreqTable.tsv: Frequency of each variant in each library. Units = per million.
* data/AAFreqTable.tsv: Frequency of each amino acid at each residue in each library. 

### ANALYSIS PIPELINE
* script/C05\_Read2Count.py
  * Input file: fastq/Wu-1\_S\*\_L001\_R1\_001.fastq
  * Output file: data/\*.count
* script/C05\_Count2VarFreq.py
  * Input file: data/\*.count
  * Output file: data/VariantFreqTable.tsv
* script/C05\_VarFreq2AAFreq.py
  * Input file: data/VariantFreqTable.tsv
  * Output file: data/AAFreqTable.tsv

### PLOTTING
* script/C05\_plot\_aaheatmap.R: Plot the frequency of each amino acid at each residue in each library as heatmaps
  * Input file: data/AAFreqTable.tsv
  * Output file: graph/HM\*.png
* script/C05\_plot\_top10varfreq.R: Plot the frequency of 10 top variants with the highest frequency in round 3 against each HA
  * Input file: data/AAFreqTable.tsv
  * Output file: graph/H?Top10Freq.png
* script/C05\_plot\_WTfrac.R: Plot the frequency of WT in each library 
  * Input file: data/AAFreqTable.tsv
  * Output file: graph/WTFrac.png
* script/C05\_plot\_H1vsH3.R: Compare the frequency of each variant in H1 round 3 vs H3 round 3
  * Input file: data/AAFreqTable.tsv
  * Output file: graph/H1H3specialists.png
