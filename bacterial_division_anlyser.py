__author__ = "Bruno Fosso"
__version__ = 2.0
__mail__ = "b.fosso@ibbe.cnr.it"

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

try:
    from Bio import SeqIO
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
except:
    sys.exit( "biopython module not found" )


def usage():
    print ('This script works in two different steps:\n'
           '\t(1) it maps the sequencing data against the Bacterial division;\n'
           '\t(2) it extract the mapped pairs and tries to map them on the human reference;\n'
           'Options:\n'
           '\t-i\tpaired-end file list: a line containing the R1, the R2 files [MANDATORY]\n'
           '\t-p\tNumber of threads/processors available [MANDATORY]'
           '\t-r\treference path [MANDATORY]\n'
           '\t-s\tscript path path [MANDATORY]\n'
           '\t-h\tprint this help.\n'
           'Usage:\n'
           '\tpython division_analyser -i read_list\n'
           '\n'
           )


script_path = ""
reference_path = ""
f_in = ""
processor = 0
try:
    opts, args = getopt.getopt( sys.argv[1:], "hi:p:r:s:" )
except getopt.GetoptError, err:
    print str( err )
    usage( )
    sys.exit( )

if len( opts ) != 0:
    for o, a in opts:
        if o == "-h":
            usage( )
            exit( )
        elif o == "-i":
            f_in = a
        elif o == "-p":
            processor = int( a )
        elif o == "-r":
            reference_path = a
        elif o == "-s":
            script_path = a
        else:
            usage( )
            sys.exit( )
else:
    usage( )
    sys.exit( )


def mapping_threads(available_processors):
    if available_processors == 0:
        print "The number of available threads/processors was not indicated!!!"
        usage( )
        exit( )
    else:
        t = 1
        if available_processors == 1 or available_processors <= 10:
            print "Due to number of available threads, single processor mapping will be performed"
            t = 1
        elif available_processors >= 15:
            t = int( processor / 15 )
        return str( t )


threads = mapping_threads( processor )


def input_data_evaluation(n, c):
    # controlliamo che il file r1 esista
    input_list = map( strip, n.split( "," ) )
    if len( input_list ) > 1:
        for r1_file in input_list:
            if os.path.exists( r1_file ) is not True:
                print "the indicated %s fastq file doesn't exist" % r1_file
                usage( )
                sys.exit( )
            else:
                try:
                    count = []
                    with open( r1_file ) as handle:
                        for title in FastqGeneralIterator( handle ):
                            count.append( title )
                            if len( count ) >= 1:
                                break
                except:
                    sys.exit( "input files are not in fastq format" )
    else:
        if os.path.exists( n ) is not True:
            print "the indicated %s fastq file doesn't exist" % n
            usage( )
            sys.exit( )
    input_list = map( strip, c.split( "," ) )
    if len( input_list ) > 1:
        for r2_file in input_list:
            if os.path.exists( r2_file ) is not True:
                print "the indicated %s fastq file doesn't exist" % r2_file
                usage( )
                sys.exit( )
            else:
                try:
                    count = []
                    with open( r2_file ) as handle:
                        for title in FastqGeneralIterator( handle ):
                            count.append( title )
                            if len( count ) >= 1:
                                break
                except:
                    sys.exit( "input files are not in fastq format" )
    else:
        if os.path.exists( c ) is not True:
            print "the indicated %s fastq file doesn't exist" % c
            usage( )
            sys.exit( )


def controllo_read_file(l):
    import os
    if os.path.exists( l ) and os.stat( l )[6] != 0:
        input_file = open( l ).readlines( )[0].split( "\t" )
        if len( input_file ) == 2:
            print "correctly formatted read_list file"
        else:
            print "uncorrectly formatted read_list file"
            exit( )
    else:
        print "read-file doesn't exist or is an empty file"
        usage( )
        exit( )


def error_file_check(l):
    from string import find
    count = 0
    for linea in open( l ):
        linea = linea.strip( )
        if find( linea.lower( ), "error" ) != -1:
            count += 1
    return count


def controll_mapping_procedure(folder):
    if os.path.exists( os.path.join( folder, "result_file.lst" ) ) and \
                    os.stat( os.path.join( folder, "result_file.lst" ) )[6] != 0:
        return os.path.join( folder, "result_file.lst" )
    else:
        return 0

def key_name(acc):
    if len(acc.split("/")) == 2:
        parts = acc.split("/")
    elif len(acc.split()) == 2:
        parts = acc.split()
    else:
        parts = acc.split()
    seq_id = parts[0]
    return seq_id

def cigar_parsing(q_name, q_len, cigar_data, value, r_name):
    mm = 0
    i = 0
    d = 0
    for item in cigar_data:
        if item[0] == 0:
            mm += item[1]
        elif item[0] == 1:
            i += item[1]
        elif item[0] == 2:
            d += item[1]
    align_len = mm + i + d
    query_aligned_len = float( mm + i )
    if query_aligned_len / q_len >= 0.7:
        if align_len != 0 and value >= 0:
            identity_percentage = ((align_len - value) / align_len) * 100
            if identity_percentage >= 93:
                evaluated_data = [q_name, r_name, fpformat.fix( identity_percentage, 2 ), str( value )]
                if align.is_read1:
                    evaluated_data.append( "r1\n" )
                elif align.is_read2:
                    evaluated_data.append( "r2\n" )
                return evaluated_data
            else:
                return None

if script_path == "":
    print "The parameter file lacks of the script_path info"
    sys.exit( )
if reference_path == "":
    print "The parameter file lacks of the reference_path info"
    sys.exit( )

# multiple mapping procedure
split2index = {}
with open( os.path.join(reference_path, "Bacteria/bacterial_bowtie_index_reference.lst") ) as a:
    for line in a:
        s = map( strip, line.split( ) )
        split2index[s[0]] = os.path.join(reference_path, s[1])

R1 = []
R2 = []
for line in open( f_in ):
    s = map( strip, line.split( "\t" ) )
    R1.append( s[0] )
    R2.append( s[1] )

if len( R1 ) != len( R2 ):
    print "The read list file contains a number of R1 files different to R2 one"
    usage( )

R1_string = ""
R2_string = ""
if len( R1 ) > 1:
    R1_string = ",".join( R1 )
    R2_string = ",".join( R2 )
else:
    R1_string = R1[0]
    R2_string = R2[0]

input_data_evaluation( R1_string, R2_string )

bacterial_folder = os.path.join( os.getcwd( ), "bacteria" )
if os.path.exists( bacterial_folder ) is False:
    os.mkdir( bacterial_folder )

split2result = {}
# mapping step
mapping_process_pid = {}
for split in split2index.keys( ):
    mapping_folder = os.path.join( bacterial_folder, split )
    if os.path.exists( mapping_folder ) is False:
        os.mkdir( mapping_folder )
        database = split2index[split]
        cmd = shlex.split( "python %s -1 %s -2 %s -d %s -s %s -f %s -t %s" % (os.path.join( script_path, "bacterial_mapper.py" ), R1_string, R2_string, database, split, mapping_folder, threads) )
        p = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
        mapping_process_pid.setdefault( split, [] )
        mapping_process_pid[split].append( p.pid )
        mapping_process_pid[split].append( mapping_folder )
    else:
        result = controll_mapping_procedure( mapping_folder )
        if result == 0:
            database = split2index[split]
            cmd = shlex.split( "python %s -1 %s -2 %s -d %s -s %s -f %s -t %s" % (os.path.join( script_path, "bacterial_mapper.py" ), R1_string, R2_string, database, split, mapping_folder, threads) )
            p = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
            mapping_process_pid.setdefault( split, [] )
            mapping_process_pid[split].append( p.pid )
            mapping_process_pid[split].append( mapping_folder )
        else:
            mapping_process_pid.setdefault( split, [] )
            mapping_process_pid[split].append( -1 )
            mapping_process_pid[split].append( mapping_folder )

completed = set( )
pe_data = os.path.join( bacterial_folder, "PE_data" )
if os.path.exists( pe_data ) is not True:
    os.mkdir( pe_data )
query2file_name = {}
counter = 1
while len( completed ) != len( mapping_process_pid ):
    for split in mapping_process_pid.keys( ):
        # print split
        pid = mapping_process_pid[split][0]
        mapping_folder = mapping_process_pid[split][1]
        process_status = ""
        if psutil.pid_exists( pid ):
            process_status = psutil.Process( pid ).status( )
        # print process_status
        if psutil.pid_exists( pid ) is False or process_status.lower( ) in ["finished", "zombie"]:
            if split not in completed:
                result = controll_mapping_procedure( mapping_folder )
                if result == 0:
                    print "%s bacterial section require a new mapping procedure" % split
                    database = split2index[split]
                    cmd = shlex.split(
                        "python %s -1 %s -2 %s -d %s -s %s -f %s -t %s" % (os.path.join( script_path, "bacterial_mapper.py" ), R1_string, R2_string, database, split, mapping_folder, threads) )
                    p = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
                    del mapping_process_pid[split]
                    mapping_process_pid.setdefault( split, [] )
                    mapping_process_pid[split].append( p.pid )
                    mapping_process_pid[split].append( mapping_folder )
                else:
                    file_name = open( result ).readlines( )[0].strip( )
                    if file_name.endswith( "bam" ):
                        sam = Samfile( file_name, "rb" )
                    elif file_name.endswith( "sam" ):
                        sam = Samfile( file_name )
                    for align in sam.fetch( until_eof=True ):
                        if align.tid != -1:
                            query_name = align.qname  # accession della read
                            query_len = float( align.rlen )  # lunghezza delle read in esame
                            # query_aligned_len = float(align.qlen)  # porzione della read allineata
                            ref_name = sam.getrname( align.tid )
                            # print query_name, ref_name
                            # print align.cigar
                            cigar = list( align.cigar )
                            nm = -1
                            for coppia in align.tags:
                                if coppia[0] == "NM":
                                    nm = float( coppia[1] )
                            stringa = cigar_parsing( query_name, query_len, cigar, nm, ref_name )
                            if stringa is not None:
                                if query2file_name.has_key( query_name ):
                                    name = os.path.join( pe_data, query2file_name[query_name] )
                                    with open( name, "a" ) as a:
                                        a.write( "\t".join( stringa ) )
                                else:
                                    query2file_name[query_name] = str( counter )
                                    counter += 1
                                    name = os.path.join( pe_data, query2file_name[query_name] )
                                    a = open( name, "w" )
                                    a.write( "\t".join( stringa ) )
                                    a.close( )
                    sam.close( )
                    completed.add( split )

match_file = open( os.path.join( bacterial_folder, "mapped_on_bacteria_total.txt" ), "w" )
print "splitted files analysis"
for name in os.listdir( pe_data ):
    # print name
    name = os.path.join( pe_data, name )
    r1_match = {}
    r2_match = {}
    with open( name ) as a:
        for line in a:
            s = map( strip, line.split( "\t" ) )
            query_name, ref_name, paired_perc_id, nm, read = s[0], s[1], float( s[2] ), float( s[3] ), s[4]
            if read == "r1":
                r1_match.setdefault( ref_name, [] )
                r1_match[ref_name].append( paired_perc_id )
                #r1_match[ref_name].append( nm )
            elif read == "r2":
                r2_match.setdefault( ref_name, [] )
                r2_match[ref_name].append( paired_perc_id )
                #r2_match[ref_name].append( nm )
    r1_acc = set( r1_match.keys( ) )
    r2_acc = set( r2_match.keys( ) )
    common = r1_acc.intersection( r2_acc )
    print common

    if len( common ) > 0:
        division_match = {}
        for ref in common:
                division_match.setdefault( ref, [] )
                perc_id = numpy.mean( [r1_match[ref][0], r2_match[ref][0]] )
                # nm_comb = numpy.mean( [r1_match[ref][1], r2_match[ref][1]] )
                division_match[ref].append( perc_id )
                # division_match[ref_name].append( nm_comb )

        # questa versione fa una riduzione dei match da prendere in considerazione sulla base degli nm
        # va modificata la struttura del dizionario!!!
        # if len( division_match.keys( ) ) != 0:
        #     for query in division_match.keys( ):
        #         over_97 = {}
        #         nm_over_97 = []
        #         for ref in division_match[query].keys( ):
        #             if division_match[query][ref][0] >= 97:
        #                 over_97[ref] = division_match[query][ref][1]
        #                 nm_over_97.append( division_match[query][ref][1] )
        #         match_list = set( )
        #         if len( nm_over_97 ) > 0:
        #             threshold = min( nm_over_97 ) + 1
        #             for ref in over_97.keys( ):
        #                 if over_97[ref] <= threshold:
        #                     match_list.add( ref )
        #         if len( match_list ) != 0:
        #             print >> match_file, query + " " + " ".join( match_list )

        if len( division_match.keys( ) ) != 0:
            match_list = set()
            match_file.add(query_name)
            for ref in division_match.keys( ):
                if division_match[ref][0] >= 97:
                    match_list.add( ref )
            if len( match_list ) > 1:
                print >> match_file, " ".join( match_list )
match_file.close( )
