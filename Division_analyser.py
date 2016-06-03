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
import time

try:
    import numpy
except:
    sys.exit( "numpy module not found" )
try:
    from pysam import *
except:
    sys.exit( "pysam module not found" )
try:
    from Bio import SeqIO
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
except:
    sys.exit( "biopython module not found" )


def usage():
    print ('This script works in three different steps:\n'
           '\t(1) it maps the sequencing data against the selected reference division;\n'
           '\t(2) it filters the mapping data\n'
           'Options:\n'
           '\t-i\tpaired-end file list: a line containing the R1, the R2 files tab separeted [MANDATORY]\n'
           '\t-d\tDivision of interest. Allowed division are: virus, bacteria, fungi, protist [MANDATORY]\n'
           '\t-r\treference path [MANDATORY]\n'
           '\t-h\tprint this help.\n'
           'Usage:\n'
           '\tpython division_analyser -i read_list -d division\n'
           '\n'
           )


try:
    opts, args = getopt.getopt( sys.argv[1:], "hi:d:" )
except getopt.GetoptError, err:
    print str( err )
    usage( )
    sys.exit( )

f_in = ""
division = ""
reference_path = ""
if len( opts ) != 0:
    for o, a in opts:
        if o == "-h":
            usage( )
            sys.exit( )
        elif o == "-i":
            f_in = a
        elif o == "-d":
            division = a.lower( )
        elif o == "-r":
            reference_path = a
        else:
            usage( )
            sys.exit( "Unhandeled option" )
else:
    usage( )
    sys.exit( )


def option_control(l, m, path ):
    if l == "":
        sys.exit( "no read list file" )
    if os.path.exists( l ) is False:
        sys.exit( "The indicated read list file doesn't exists" )
    if m == "":
        sys.exit( "Missing Division choice" )
    if os.path.exists( os.path.dirname( m ) ) is False:
        sys.exit( "The indicated division data don't exists" )
    if path == "":
        sys.exit( "Missing refence path" )


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


division2bowtie_index_files = {"virus": os.path.join(reference_path, "Virus/bowtie_index/Viral_division"),
                               "fungi": os.path.join(reference_path, "Fungi/bowtie_index/Fungal_division"),
                               "protist": os.path.join(reference_path, "Protist/bowtie_index/Protista_division")}

option_control( f_in, division2bowtie_index_files[division], reference_path )

script_time = time.strftime( "%d/%m/%Y %H:%M:%S", time.localtime( time.time( ) ) )
print "START", script_time

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


def key_name(acc):
    if len( acc.split( "/" ) ) == 2:
        parts = acc.split( "/" )
    elif len( acc.split( ) ) == 2:
        parts = acc.split( )
    else:
        parts = acc.split( )
    seq_id = parts[0]
    return seq_id


def error_file_check(l):
    from string import find
    count = 0
    for line_dat in open( l ):
        line_dat = line_dat.strip( )
        if find( line_dat.lower( ), "error" ) != -1:
            count += 1
    return count


input_data_evaluation( R1_string, R2_string )

wd = os.getcwd( )
if os.path.exists( os.path.join( wd, division ) ):
    pass
else:
    os.mkdir( os.path.join( wd, division ) )

# esecuzione del mapping divisione specifico
print "bowtie on division data"
sam_output = os.path.join( wd, division, "mapping_on_%s.sam" % division )
database = division2bowtie_index_files[division]
cmd = shlex.split(
    "bowtie2 -q -N1 -k100 -L20 --no-unal -x  %s -1 %s -2 %s -S %s -p 20" % (database, R1_string, R2_string, sam_output) )
tmp = open( os.path.join( wd, division, "error.lst" ), "w" )
p = subprocess.Popen( cmd, stderr=tmp )
p.wait( )
tmp.close( )
# sam conversion to bam
if error_file_check( os.path.join( wd, division, "error.lst" ) ) == 0:
    bam_output = os.path.join( wd, division, "mapping_on_%s.bam" % division )
    cmd = shlex.split( "samtools view -bS %s -o %s" % (sam_output, bam_output) )
    tmp = open( os.path.join( wd, division, "error.lst" ), "w" )
    p = subprocess.Popen( cmd, stderr=tmp )
    p.wait( )
    tmp.close( )
else:
    print "Errors during the mapping procedure for the %s division" % division
    print "please see error.lst file"
    sys.exit( )

if error_file_check( os.path.join( wd, division, "error.lst" ) ) == 0:
    sam = Samfile( bam_output, "rb" )
    # os.remove(sam_output)
else:
    print "Errors happened during the samtools view procedure"
    print "BAM conversion was stopped. Please control the error.lst file"
    print "SAM will be used to identity significant mapping data. This could require long processing time"
    sam = Samfile( sam_output )
division_match = {}
# parsing dei dati ottenuti
print "sam parsing"
script_time = time.strftime( "%d/%m/%Y %H:%M:%S", time.localtime( time.time( ) ) )
print "START", script_time
pe_data = os.path.join( wd, division, "PE_data" )
if os.path.exists( pe_data ) is not True:
    os.mkdir( pe_data )
query2file_name = {}
counter = 1
for align in sam:
    if align.tid != -1:
        query_name = align.qname  # accession della read
        query_len = float( align.rlen )  # lunghezza delle read in esame
        # query_aligned_len = float(align.qlen)  # porzione della read allineata
        ref_name = sam.getrname( align.tid )
        # print query_name, ref_name
        cigar = list( align.cigar )
        print query_name, query_len, ref_name, cigar
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

match_file = open( os.path.join( division, "mapped_on_" + division + "_total.txt" ), "w" )
print "splitted files analysis"
for name in os.listdir( pe_data ):
    # print name
    name = os.path.join( pe_data, name )
    r1_match = {}
    r2_match = {}
    division_match = {}
    with open( name ) as a:
        for line in a:
            s = map( strip, line.split( "\t" ) )
            query_name, ref_name, paired_perc_id, nm, read = s[0], s[1], float( s[2] ), float( s[3] ), s[4]
            if read == "r1":
                r1_match.setdefault( query_name, {} )
                r1_match[query_name].setdefault( ref_name, [] )
                r1_match[query_name][ref_name].append( paired_perc_id )
                r1_match[query_name][ref_name].append( nm )
            elif read == "r2":
                r2_match.setdefault( query_name, {} )
                r2_match[query_name].setdefault( ref_name, [] )
                r2_match[query_name][ref_name].append( paired_perc_id )
                r2_match[query_name][ref_name].append( nm )
    r1_acc = set( r1_match.keys( ) )
    r2_acc = set( r2_match.keys( ) )
    common = r1_acc.intersection( r2_acc )

    for key in common:
        # definiamo le coppie per cui le due read mappano sulla stessa seq di riferimento
        common_match = set( r1_match[key].keys( ) ).intersection( set( r2_match[key].keys( ) ) )
        if len( common_match ) > 0:
            for ref_name in common_match:
                division_match.setdefault( key, {} )
                division_match[key].setdefault( ref_name, [] )
                perc_id = numpy.mean( [r1_match[key][ref_name][0], r2_match[key][ref_name][0]] )
                nm_comb = numpy.mean( [r1_match[key][ref_name][1], r2_match[key][ref_name][1]] )
                division_match[key][ref_name].append( perc_id )
                division_match[key][ref_name].append( nm_comb )

    # questa versione fa una riduzione dei match da prendere in considerazione sulla base degli nm
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
        for query in division_match.keys( ):
            match_list = set( )
            for ref in division_match[query].keys( ):
                if division_match[query][ref][0] >= 97:
                    match_list.add( ref )
            if len( match_list ) != 0:
                print >> match_file, query + " " + " ".join( match_list )
match_file.close( )

script_time = time.strftime( "%d/%m/%Y %H:%M:%S", time.localtime( time.time( ) ) )
print "END", script_time
