__author__ = "Bruno Fosso"
__version__ = 2.0
__mail__ = "b.fosso@ibbe.cnr.it"

import fpformat
import getopt
import os
import shlex
import subprocess
import sys
import shutil
import time
from string import strip

try:
    import psutil
except:
    sys.exit( "psutil module not found" )
try:
    import numpy
except:
    sys.exit( "numpy module not found" )
try:
    from Bio import SeqIO
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
except:
    sys.exit( "biopython module not found" )


def usage():
    print ('This script extracts the PE reads belonging to a specific taxa:\n'
           '\tit requires in input a single or a list of NCBI taxonomy identifier\n'
           'Options:\n'
           '\t-t\tNCBI taxonomy ID. To extract human sequences please use 9606\n'
           '\t-l\ta text file containing a list of NCBI taxonomy ID, one per line\n'
           '\t-u\textract unassigned PE reads (the ambiguous PE reads are not considered)\n'
           '\t-a\textract ambiguous PE reads\n'
           '\t-h\tprint this help.\n'
           'Usage:\n'
           '\tpython PE_extraction.py -t 9606\n'
           '\n'
           )




taxid_list = set( )
taxid_file = ""
try:
    opts, args = getopt.getopt( sys.argv[1:], "ht:l:ua" )
except getopt.GetoptError, err:
    print str( err )
    usage( )
    sys.exit( )

if len( opts ) != 0:
    for o, a in opts:
        if o == "-h":
            usage( )
            exit( )
        elif o == "-t":
            taxid_list.add( a )
        elif o == "-u":
            taxid_list.add( "unassigned" )
        elif o == "-a":
            taxid_list.add( "ambiguous" )
        elif o == "-l":
            taxid_file = a
        else:
            usage( )
            sys.exit( )
else:
    usage( )
    sys.exit( )


def populate_taxid_list(name, lista):
    if os.path.exists( name ):
        if os.stat( name )[6] != 0:
            with open( name ) as l:
                for linea in l:
                    lista.add( linea.strip( ) )
            return lista
        else:
            print "The NCBI list file is empty"
            print usage( )
            sys.exit( )
    else:
        print "The NCBI list file doesn't exist"
        print usage( )
        sys.exit( )


def key_name(acc_num):
    if len( acc_num.split( "/" ) ) == 2:
        parts = acc_num.split( "/" )
    elif len( acc_num.split( ) ) == 2:
        parts = acc_num.split( )
    else:
        parts = acc_num.split( )
    seq_id = parts[0]
    return seq_id


if taxid_file is not "":
    taxid_list = populate_taxid_list( taxid_file, taxid_list )
else:
    if len( taxid_list ) < 1:
        print "neither the option nor b have been indicated"
        print usage( )
        sys.exit( )

wd = os.getcwd( )
# list file creation
acc2taxid = {}
total_assigned = set()
for division in ["bacteria", "virus", "protist", "fungi"]:
    print division
    result_file = os.path.join( wd, division, "%s_classification_data.txt_0.5" % division )
    if os.path.exists( result_file ) and os.stat( result_file )[6] != 0:
        with open( result_file ) as m:
            for line in m:
                field = map( strip, line.split( "\t" ) )
                acc, path = field[0], set( field[2].split( ";" ) )
                extractable = path.intersection( taxid_list )
                total_assigned.add(acc)
                if len( extractable ) >= 1:
                    acc2taxid[acc] = extractable
if "9606" in taxid_list:
    acc2taxid.setdefault( "9606", set( ) )
    with open( os.path.join( wd, "mapping_on_human", "only_human.lst" ) ) as a:
        for line in a:
            acc2taxid["9606"].add( line.strip( ) )
            total_assigned.add( line.strip( ) )
if "ambiguous" in taxid_list:
    acc2taxid.setdefault( "ambiguous", set( ) )
    with open( os.path.join( wd, "ambiguos_pe_read.lst" ) ) as a:
        for line in a:
            acc2taxid["ambiguous"].add( line.strip( ) )
            total_assigned.add( line.strip( ) )
print "A total number of %i was selected for the extraction" % len( acc2taxid )

script_time = time.strftime("%d_%m_%Y_%H_%M_%S", time.localtime(time.time()))
extraction_folder = os.path.join( wd, "taxid_PE_file_" + script_time )
if os.path.exists( extraction_folder ):
    shutil.move( extraction_folder, extraction_folder + "_old" )
    os.mkdir( extraction_folder )
else:
    os.mkdir( extraction_folder )

wrote = 0
taxid2seq_data = {}
with open("read_list_cleaned") as file_list:
    for line in file_list:
        s = map( strip, line.split( "\t" ) )
        R1 = s[0]
        R2 = s[1]
        with open( R1 ) as fastq:
            for title, seq, qual in FastqGeneralIterator( fastq ):
                stringa = "@%s\n%s\n+\n%s\n" % (title, seq, qual)
                if acc2taxid.has_key( key_name( title ) ):
                    for taxid in acc2taxid[key_name( title )]:
                        seq_name = os.path.join( extraction_folder, "%s_R1.fastq" % taxid )
                        if os.path.exists( seq_name ):
                            with open( name, "a" ) as a:
                                a.write( stringa )
                        else:
                            a = open( seq_name, "w" )
                            a.write( stringa )
                            a.close( )
                else:
                    if "unassigned" in taxid_list:
                        seq_name = os.path.join( extraction_folder, "unassigned_R1.fastq" )
                        if os.path.exists( seq_name ):
                            with open( seq_name, "a" ) as a:
                                a.write( stringa )
                        else:
                            a = open( seq_name, "w" )
                            a.write( stringa )
                            a.close( )
        with open( R2 ) as fastq:
            for title, seq, qual in FastqGeneralIterator( fastq ):
                stringa = "@%s\n%s\n+\n%s\n" % (title, seq, qual)
                if acc2taxid.has_key( key_name( title ) ):
                    for taxid in acc2taxid[key_name( title )]:
                        seq_name = os.path.join( extraction_folder, "%s_R2.fastq" % taxid )
                        if os.path.exists( seq_name ):
                            with open( seq_name, "a" ) as a:
                                a.write( stringa )
                        else:
                            a = open( seq_name, "w" )
                            a.write( stringa )
                            a.close( )
                else:
                    if "unassigned" in taxid_list:
                        seq_name = os.path.join( extraction_folder, "unassigned_R2.fastq" )
                        if os.path.exists( seq_name ):
                            with open( seq_name, "a" ) as a:
                                a.write( stringa )
                        else:
                            a = open( seq_name, "w" )
                            a.write( stringa )
                            a.close( )
