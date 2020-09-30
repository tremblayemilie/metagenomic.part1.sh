# metagenomic.part1.sh
Ion Torrent Raw data process until swithching to part 1 in R

barcodes_Debbie_pollen.txt is a text file with tab separated columns including data associated with your samples. 
Example:

#SampleID       Barcode ID      BarcodeSequence Primer  CollectionDate
DPS01   ITS1F_BC01      CTAAGGTAAC      GATCTTGGTCATTTAGAGGAAGTAA       ITS1F   2017-05-17
DPS02   ITS1F_BC02      TAAGGAGAAC      GATCTTGGTCATTTAGAGGAAGTAA       ITS1F   2017-05-27

#the script fastqtofasta.py can be found online (https://github.com/chanzuckerberg/shasta/blob/master/scripts/FastqToFasta.py)

or here:
#!/usr/bin/python3

import os
import sys

helpMessage = """
This script will convert a fastq file to fasta.
Invoke with two arguments:
- Name of the input fastq file.
- Name of the output fasta file.
"""

if not len(sys.argv)==3:
    print(helpMessage)
    exit(1) 
    
inputFileName = sys.argv[1]
outputFileName = sys.argv[2]

command = 'cat ' + inputFileName + ' | awk \'{if(NR%4==1) {printf(">%s\\n",substr($0,2));} else if(NR%4==2) print;}\' > ' + outputFileName
    
os.system(command)    
