```
usage: fastq_preprocessor.py [-h] [--timestamp] [--newDIR NEWDIR]
                             [--newDir NEWDIR] [--moveRaw MOVERAW]
                             [--rename RENAME] [--force FORCE]
                             [--checkMatch CHECKMATCH] [--debug]
                             [--NCORE NCORE]
                             samplePATH

# Usage: (python) preprocessor.py /path/to/FASTQ_DIR
# Example: preprocessor.py /media/pw_synology3/PW_HiSeq_data/RNA-seq/Raw_data/testONLY/133R/BdPIFs-32747730/133E_23_DN-40235206
# Purpose: Download fastq files from the supplied path and 
#    combine them into R1.fastq and R2.fastq
# 
# Created:7  OCT 2016, Hui@SLCU map-RNA-seq.py
# Update: 11 Feb 2017. Hui@SLCU
# Update: 29 May 2018. Feng@SLCU preprocessor.py
# Update: 08 Apr 2020. Feng

Example Usage:
    fastq_preprocess test_data --newDIR test_out
    preprocessor.py test_data --newDIR test_out
Install:
    python -m pip install "pip>=19.0" --upgrade
    python -m pip install fastq_preprocessor@https://github.com/shouldsee/fastq_preprocessor/tarball/0.0.2

CHANGELOG:
# 0.0.2
- fixed a bug for concatenating fastq files

positional arguments:
  samplePATH

optional arguments:
  -h, --help            show this help message and exit
  --timestamp
  --newDIR NEWDIR
  --newDir NEWDIR
  --moveRaw MOVERAW
  --rename RENAME
  --force FORCE
  --checkMatch CHECKMATCH
  --debug
  --NCORE NCORE
```
