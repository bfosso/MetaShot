__author__ = "Bruno Fosso"
__version__ = 2.0

import getopt
import gzip
import os
import shlex
import subprocess
import sys
import time

try:
    import numpy
except:
    sys.exit("numpy module not found")
try:
    from pysam import Samfile
except:
    sys.exit("pysam module not found")
try:
    from Bio import SeqIO
except:
    sys.exit("biopython module not found")
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def usage():
    print ('This script works finds and remove phage PhiX sequences:\n'
           'Options:\n'
           '\t-1\tR1.fastq file [MANDATORY]\n'
           '\t-2\tR2.fastq file [MANDATORY]\n'
           '\t-o\toutput folder [MANDATORY]\n'
           '\t-p\treference path [MANDATORY]\n'
           '\t-h\tprint this help.\n'
           'Usage:\n'
           '\tpython Phix_cleaner.py -1 R1.fastq -2 R2.fastq -o output_folder\n'
           '\n'
           )


try:
    opts, args = getopt.getopt(sys.argv[1:], "h1:2:o:p:")
except getopt.GetoptError, err:
    print str(err)
    usage()
    sys.exit()

R1 = ""
R2 = ""
output_folder = ""
reference_path = ""
if len(opts) != 0:
    for o, a in opts:
        if o == "-h":
            usage()
            sys.exit()
        elif o == "-1":
            R1 = a
        elif o == "-2":
            R2 = a
        elif o == "-o":
            output_folder = a
        elif o == "-p":
            reference_path = a
        else:
            usage()
            sys.exit("Unhandeled option")
else:
    usage()
    sys.exit()


def control_options(aa, b, folder):
    if aa == "":
        sys.exit("no R1 file")
    if b == "":
        sys.exit("no R2 file")
    if folder == "":
        sys.exit("no output folder")


def return_acc(acc):
    parts = acc.split()[0].split("/")
    return parts[0]


script_time = time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
print "START", script_time

control_options(R1, R2, output_folder)

phix_folder = os.path.join(output_folder)
if os.path.exists(phix_folder):
    pass
else:
    os.mkdir(phix_folder)

database = os.path.join(reference_path, "Phix/phix")
sam_output = os.path.join(phix_folder, "phix_mapping_data.sam")
cmd = shlex.split("bowtie2 -q --no-unal -x  %s -1 %s -2 %s -S %s -p 5 --mm" % (database, R1, R2, sam_output))
# p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
p = subprocess.Popen(cmd)
p.wait()

sam = Samfile(sam_output)
phix_reads = set()
for align in sam.fetch(until_eof=True):
    if align.tid != -1:
        query_name = align.qname  # accession della read
        query_len = float(align.rlen)  # lunghezza delle read in esame
        threshold = int((query_len / 100) * 3) + 1
        if align.cigar is not None:
            nm = -1
            for coppia in align.tags:
                if coppia[0] == "NM":
                    nm = float(coppia[1])
            if 0 <= nm <= threshold:
                phix_reads.add(return_acc(query_name))

print len(phix_reads)
result = open(os.path.join(phix_folder, "Phix_sequences.lst"), "w")
result.write("\n".join(list(phix_reads)))
result.close()

#################
# DATA CLEANING #
#################
cleaned_r1 = open(os.path.join(phix_folder, "R1_no_phix.fastq"), "w")
if R1.endswith("gz"):
    l1 = gzip.open(R1)
else:
    l1 = open(R1)
for title, seq, qual in FastqGeneralIterator(l1):
    if return_acc(title) not in phix_reads:
        cleaned_r1.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
l1.close()
cleaned_r1.close()

cleaned_r2 = open(os.path.join(phix_folder, "R2_no_phix.fastq"), "w")
if R1.endswith("gz"):
    l2 = gzip.open(R2)
else:
    l2 = open(R2)
for title, seq, qual in FastqGeneralIterator(l2):
    if return_acc(title) not in phix_reads:
        cleaned_r2.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
l2.close()
cleaned_r2.close()

script_time = time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(time.time()))
print "END", script_time
