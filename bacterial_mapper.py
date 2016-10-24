__author__ = "Bruno Fosso"
__version__ = 2.0
import fpformat
import getopt
import os
import shlex
import subprocess
import sys
from string import strip
import psutil
import numpy
from pysam import *


def usage():
    print ('This script performs the mapping on the sections of the bacterial division:\n'
           'Options:\n'
           '\t1\tfastq file containing the R1 sequences[MANDATORY]\n'
           '\t2\tfastq file containing the R2 sequences[MANDATORY]\n'
           '\td\tsection of the bacterial division[MANDATORY]\n'
           '\ts\tbacterial split[MANDATORY]\n'
           '\tf\tworking folder[MANDATORY]\n'
           '\tt\tnumber of available threads[MANDATORY]\n'
           '\t-h\tprint this help.\n'
           'Usage:\n'
           '\tpython bacterial_mapper.py -1 r1.fastq -2 r2.fastq -d bacteria_section_index -t 10 -f working folder\n'
           '\n'
           )


threads = 0
R1 = ""
R2 = ""
database = ""
mapping_folder = ""
split = ""
try:
    opts, args = getopt.getopt( sys.argv[1:], "h1:2:d:t:f:s:" )
except getopt.GetoptError, err:
    print str( err )
    usage( )
    sys.exit( )
if len( opts ) != 0:
    for o, a in opts:
        if o == "-h":
            usage( )
            sys.exit( )
        elif o == "-1":
            R1 = a
        elif o == "-2":
            R2 = a
        elif o == "-d":
            database = a
        elif o == "-t":
            threads = a
        elif o == "-f":
            mapping_folder = a
        elif o == "-s":
            split = a
        else:
            usage( )
            sys.exit( )


def option_control(l, m, n, c):
    print l, m, n, c
    # controlliamo che il file r1 esista
    input_list = map( strip, l.split( "," ) )
    if len( input_list ) > 1:
        for r1_file in input_list:
            if os.path.exists( r1_file ) is not True:
                print "the indicated %s fastq file doesn't exist" % r1_file
                usage( )
                sys.exit( )
    else:
        if os.path.exists( l ) is not True:
            print "the indicated %s fastq file doesn't exist" % l
            usage( )
            sys.exit( )
    input_list = map( strip, m.split( "," ) )
    if len( input_list ) > 1:
        for r2_file in input_list:
            if os.path.exists( r2_file ) is not True:
                print "the indicated %s fastq file doesn't exist" % r2_file
                usage( )
                sys.exit( )
    else:
        if os.path.exists( m ) is not True:
            print "the indicated %s fastq file doesn't exist" % m
            usage( )
            sys.exit( )
    if type( n ) == int or n == "0":
        usage( )
        sys.exit( )
    if mapping_folder == "":
        print "No mapping folder"
        usage( )
        sys.exit( )


def error_file_check(l):
    from string import find
    count = 0
    for line in open( l ):
        line = line.strip( )
        if find( line.lower( ), "error" ) != -1:
            count += 1
    return count


option_control( R1, R2, threads, mapping_folder )

sam_output = os.path.join( mapping_folder, "mapping_on_bacterial_%s.sam" % split )
cmd = shlex.split(
    "bowtie2 -q -N1 -k50 -L20 --no-unal -x  %s -1 %s -2 %s -S %s -p %s --mm" % (database, R1, R2, sam_output, threads) )
tmp = open( os.path.join( mapping_folder, "error.lst" ), "w" )
p = psutil.Popen( cmd, stderr=tmp )
p.wait( )
tmp.close( )
# sam conversion to bam
if error_file_check( os.path.join( mapping_folder, "error.lst" ) ) == 0:
    bam_output = os.path.join( mapping_folder, "mapping_on_bacterial_%s.bam" % split )
    cmd = shlex.split( "samtools view -bS %s -o %s" % (sam_output, bam_output) )
    tmp = open( "error.lst", "w" )
    p = psutil.Popen( cmd, stderr=tmp )
    p.wait( )
    tmp.close( )
else:
    print "Errors during the mapping procedure for the %s section of Bacterial database" % split
    print "please see error.lst file"
    sys.exit( )
if error_file_check( "error.lst" ) == 0:
    tmp = open( os.path.join( mapping_folder, "result_file.lst" ), "w" )
    tmp.write( bam_output )
    tmp.close( )
    os.remove( sam_output )
else:
    print "Errors happened during the samtools view procedure"
    print "BAM conversion was stopped. Please control the error.lst file"
    print "SAM will be used to identity significant mapping data. This could require long processing time"
    tmp = open( os.path.join( mapping_folder, "result_file.lst" ), "w" )
    tmp.write( sam_output )
    tmp.close( )

sys.exit( )
