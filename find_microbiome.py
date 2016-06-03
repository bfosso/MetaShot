__author__ = "Bruno Fosso"
__version__ = 2.0
import fpformat
import getopt
import os
import shlex
import subprocess
import sys
from string import strip
import time
import psutil

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
except:
    sys.exit( "biopython module not found" )
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def usage():
    print ('This script works in two different steps:\n'
           '\t(1) it fast maps the sequencing data against all the refereces searching for all possible microbial reads;\n'
           '\t(2) it extracts all possible microbial reads\n'
           'Options:\n'
           '\t-i\tpaired-end file list: a line containing the R1, the R2 files tab separeted [MANDATORY]\n'
           '\t-r\treference path [MANDATORY]\n'
           '\t-h\tprint this help.\n'
           'Usage:\n'
           '\tpython division_analyser -i read_list -d division\n'
           '\n'
           )


try:
    opts, args = getopt.getopt( sys.argv[1:], "hi:r:" )
except getopt.GetoptError, err:
    print str( err )
    usage( )
    sys.exit( )

f_in = ""
reference_path = ""
if len( opts ) != 0:
    for o, a in opts:
        if o == "-h":
            usage( )
            sys.exit( )
        elif o == "-i":
            f_in = a
        elif o == "-r":
            reference_path = a
        else:
            usage( )
            sys.exit( "Unhandeled option" )
else:
    usage( )
    sys.exit( )


def control_options(l, ref_path):
    if l == "":
        sys.exit( "no read list file" )
    if ref_path == "":
        sys.exit( "no reference_path" )


def controllo_read_file(l):
    # verify the existens of input files
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    if os.path.exists( l ) and os.stat( l )[6] != 0:
        # read list file exists
        with open( l ) as input_file:
            field = input_file.readline( ).split( "\t" )
        if len( field ) == 2:
            read1 = field[0]
            read2 = field[1].strip( )
            print "correctly formatted read_list file"
            print "fastq control"
            if os.path.exists( read1 ) and os.stat( read1 )[6] != 0:
                pass
            else:
                sys.exit( "R1 fastq file doesn't exist or is an empty file" )
            if os.path.exists( read2 ) and os.stat( read2 )[6] != 0:
                pass
            else:
                sys.exit( "R1 fastq file doesn't exist or is an empty file" )
            try:
                from Bio.SeqIO.QualityIO import FastqGeneralIterator
                count = []
                print read1
                with open( read1 ) as handle:
                    for seq_acc_id in FastqGeneralIterator( handle ):
                        count.append( seq_acc_id )
                        if len( count ) >= 1:
                            break
            except:
                sys.exit( "input files are not in fastq format" )
            try:
                print read2
                count = []
                with open( read2 ) as handle:
                    for seq_acc_id in FastqGeneralIterator( handle ):
                        count.append( seq_acc_id )
                        if len( count ) >= 1:
                            break
            except:
                sys.exit( "input files are not in fastq format" )
        else:
            sys.exit( "uncorrectly formatted read_list file" )
    else:
        sys.exit( "read-file doesn't exist or is an empty file" )


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


def return_acc(acc):
    parts = acc.split( )[0].split( "/" )
    return parts[0]


control_options( f_in, reference_path )
bowtie_index_files = {}
for line in open( os.path.join( reference_path, "find_microbiome_index.tsv" ) ):
    s = map( strip, line.split( ) )
    bowtie_index_files[s[0]] = os.path.join( reference_path, s[1] )

controllo_read_file( f_in )

script_time = time.strftime( "%d/%m/%Y %H:%M:%S", time.localtime( time.time( ) ) )
print "START", script_time

R1 = ""
R2 = ""
for line in open( f_in ):
    line = line.strip( )
    s = line.split( "\t" )
    R1 = s[0]
    R2 = s[1]


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


def controll_mapping_procedure(sam_file):
    if os.path.exists( os.path.join( sam_file ) ) and os.stat( os.path.join( sam_file ) )[6] != 0:
        try:
            controllo_sam = Samfile( os.path.join( sam_file ) )
            controllo_sam.close( )
            risultato = 1
        except:
            risultato = 0
    else:
        risultato = 0
        return risultato


wd = os.getcwd( )
premapping_folder = os.path.join( wd, "premapping_folder" )
if os.path.exists( premapping_folder ):
    pass
else:
    os.mkdir( premapping_folder )

index2result = {}
# mapping step
mapping_process_pid = {}
for split in bowtie_index_files.keys( ):
    mapping_folder = os.path.join( premapping_folder, split )
    sam_output = os.path.join( mapping_folder, split + ".sam" )
    if os.path.exists( mapping_folder ) is False:
        os.mkdir( mapping_folder )
        database = bowtie_index_files[split]
        cmd = shlex.split( "bowtie2 -q --no-unal -x  %s -1 %s -2 %s -S %s -p 2" % (database, R1, R2, sam_output) )
        p = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
        mapping_process_pid.setdefault( split, [] )
        mapping_process_pid[split].append( p.pid )
        mapping_process_pid[split].append( mapping_folder )
        mapping_process_pid[split].append( sam_output )
    else:
        result = controll_mapping_procedure( sam_output )
        if result == 0:
            database = bowtie_index_files[split]
            cmd = shlex.split( "bowtie2 -q --no-unal -x  %s -1 %s -2 %s -S %s -p 2" % (database, R1, R2, sam_output) )
            p = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
            mapping_process_pid.setdefault( split, [] )
            mapping_process_pid[split].append( p.pid )
            mapping_process_pid[split].append( mapping_folder )
            mapping_process_pid[split].append( sam_output )

to_analyse = set( )
completed = set( )
while len( completed ) != len( mapping_process_pid ):
    for split in mapping_process_pid.keys( ):
        # print split
        pid = mapping_process_pid[split][0]
        mapping_folder = mapping_process_pid[split][1]
        sam_output = mapping_process_pid[split][2]
        process_status = ""
        if psutil.pid_exists( pid ):
            process_status = psutil.Process( pid ).status( )
            print pid, process_status
        # print process_status
        if psutil.pid_exists( pid ) is False or process_status.lower( ) == "zombie":
            if split not in completed:
                result = controll_mapping_procedure( sam_output )
                if result == 0:
                    print pid, mapping_folder, sam_output
                    database = bowtie_index_files[split]
                    cmd = shlex.split( "bowtie2 -q --no-unal -x  %s -1 %s -2 %s -S %s -p 2" % (database, R1, R2, sam_output) )
                    p = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
                    mapping_process_pid.setdefault( split, [] )
                    mapping_process_pid[split].append( p.pid )
                    mapping_process_pid[split].append( mapping_folder )
                    mapping_process_pid[split].append( sam_output )
                else:
                    bam_output = sam_output.replace( ".sam", ".bam" )
                    cmd = shlex.split( "samtools view -bS %s -o %s" % (sam_output, bam_output) )
                    tmp = open( os.path.join( mapping_folder, "error.lst" ), "w" )
                    p = subprocess.Popen( cmd, stderr=tmp )
                    p.wait( )
                    tmp.close( )
                    if error_file_check( os.path.join( mapping_folder, "error.lst" ) ) == 0:
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
                    file_name = open( os.path.join( mapping_folder, "result_file.lst" ) ).readlines( )[0].strip( )
                    if file_name.endswith( "bam" ):
                        sam = Samfile( file_name, "rb" )
                        for align in sam.fetch( until_eof=True ):
                            if align.tid != -1:
                                query_name = align.qname  # accession della read
                                query_len = float( align.rlen )  # lunghezza delle read in esame
                                threshold = int( (query_len / 100) * 9 ) + 1
                                if align.cigar is not None:
                                    nm = -1
                                    for coppia in align.tags:
                                        if coppia[0] == "NM":
                                            nm = float( coppia[1] )
                                    if 0 <= nm <= threshold:
                                        to_analyse.add( key_name( query_name ) )
                    elif file_name.endswith( "sam" ):
                        sam = Samfile( file_name )
                        for align in sam.fetch( until_eof=True ):
                            if align.tid != -1:
                                query_name = align.qname  # accession della read
                                query_len = float( align.rlen )  # lunghezza delle read in esame
                                threshold = int( (query_len / 100) * 9 ) + 1
                                if align.cigar is not None:
                                    nm = -1
                                    for coppia in align.tags:
                                        if coppia[0] == "NM":
                                            nm = float( coppia[1] )
                                    if 0 <= nm <= threshold:
                                        to_analyse.add( key_name( query_name ) )
                    completed.add( split )

print "PE da mappare", len( to_analyse )
candidate_r1 = open( os.path.join( wd, "R1_micro_candidates.fastq" ), "w" )
with open( R1 ) as fastq:
    for title, seq, qual in FastqGeneralIterator( fastq ):
        if key_name( title ) in to_analyse:
            candidate_r1.write( "@%s\n%s\n+\n%s\n" % (title, seq, qual) )
candidate_r1.close( )

candidate_r2 = open( os.path.join( wd, "R2_micro_candidates.fastq" ), "w" )
with open( R2 ) as fastq:
    for title, seq, qual in FastqGeneralIterator( fastq ):
        if key_name( title ) in to_analyse:
            candidate_r2.write( "@%s\n%s\n+\n%s\n" % (title, seq, qual) )
candidate_r2.close( )
script_time = time.strftime( "%d/%m/%Y %H:%M:%S", time.localtime( time.time( ) ) )
print "END", script_time
