__author__ = 'Bruno Fosso'
__version__ = 1.0
__manteiner__ = "Bruno Fosso"
__mail__ = "b.fosso@ibbe.cnr.it"

import getopt
import gzip
import os
import psutil
import shlex
import subprocess
import sys
from string import strip

try:
    from pysam import *
except:
    print "pysam is not installed"
    sys.exit( )
try:
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
except:
    print "biopython is not installed"
    sys.exit( )

script_info = {"Description": """
This script performs the whole MetaShot computation
MetaShot Computation can be distinguished in two principal steps:
1) Pre-processing procedures: by applying FaQCs, it removes low-quality, low complexity and short reads (under 50 nt)
2) Taxonomic investigation
""", "usage": """
MetaShot is designed to analyse Illumina Paired-End (PE) data.
\n
Script Options:
\t-m\tA textual file containing the R1 and R2 PE reads file names, tab separated [MANDATORY].
\t  \tIf a sample has been splitted in more flowcell LINES, please insert in the file a line per PE reads
\t  \tExample: The sample1 has been sequenced in 3 flowcell lines. The read file content will be the following:
\t  \tsample1_L001_R1_1.fastq\tsample1_L001_R2_1.fastq
\t  \tsample1_L002_R1_2.fastq\tsample1_L002_R2_2.fastq
\t  \tsample1_L003_R1_3.fastq\tsample1_L003_R2_3.fastq
\t-p\tparameters files: a file containing all the information required for the MetaShot application [MANDATORY]
\t-t\tNumber of available CPU [MANDATORY]
\t-x\tPhage PhiX removal. This option enable the procedure of Phix removal implemented in FaQS. It may be very slow (approximately 7-time slower)
\t-h\tprint this help
\n
Example:
\tMetaShot_Master_script.py -p param_file.txt -m read_files -t 30 -s DNA
"""}

script_path = "/home/bfosso/share/MetaShot_reference/script/MetaShot_execution"

multiple_input_data = []
processor_number = ""
wd = os.getcwd( )
phix = ""
parameters_file = ""
try:
    opts, args = getopt.getopt( sys.argv[1:], "m:t:p:xh" )
except getopt.GetoptError, err:
    print str( err )
    print script_info["Description"]
    sys.exit( )
if len( opts ) != 0:
    for o, a in opts:
        if o == "-m":
            with open( a ) as read_file:
                for line in read_file:
                    s = map( strip, line.split( ) )
                    print s
                    if len( s ) == 2:
                        multiple_input_data.append( s )
            if len( multiple_input_data ) == 0:
                print "ERROR: not correctly formatted input file"
                print script_info["usage"]
                sys.exit( )
        elif o == "-t":
            processor_number = a
        elif o == "-x":
            phix = "-phiX"
        elif o == "-p":
            parameters_file = a
        elif o == "-h":
            print script_info["usage"]
            sys.exit( )
        else:
            print script_info["usage"]
            sys.exit( )
else:
    print script_info["usage"]
    sys.exit( )


#######################
# function definition #
#######################
def pid_status(process_pid):
    if psutil.pid_exists( process_pid ):
        status = psutil.Process( process_pid ).status( )
    else:
        status = "finished"
    return status


def verify_input_parameters(n_pr, p_file):
    if n_pr == "":
        print "Error!!! The -t option has not been inserted."
        sys.exit( )
    if p_file == "":
        print "Error!!! The -p option has not been inserted."
        print script_info["usage"]
        sys.exit( )
    else:
        if os.path.exists( p_file ):
            if os.stat( p_file )[6] == 0:
                print "Error!!! The indicated parameter file is empty"
                print script_info["usage"]
                sys.exit( )
        else:
            print "Error. The selected parameter file doesn't exist"
            print script_info["usage"]
            sys.exit( )


def control_mapping_procedure(sam_file):
    if os.path.exists( os.path.join( sam_file ) ) and os.stat( os.path.join( sam_file ) )[6] != 0:
        try:
            controllo_sam = Samfile( os.path.join( sam_file ) )
            controllo_sam.close( )
            result_number = 1
        except:
            result_number = 0
    else:
        result_number = 0
    return result_number


def verify_fastq(fastq1, fastq2):
    count_1 = 0
    if fastq1.endswith( "gz" ):
        l1 = gzip.open( fastq1 )
    else:
        l1 = open( fastq1 )
    for title, seq, qual in FastqGeneralIterator( l1 ):
        count_1 += 1
    l1.close( )
    count_2 = 0
    if fastq2.endswith( "gz" ):
        l2 = gzip.open( fastq2 )
    else:
        l2 = open( fastq2 )
    for acc, sequence, quality in FastqGeneralIterator( l2 ):
        count_2 += 1
    l2.close( )
    if count_1 == count_2:
        fastq_result = 1
    else:
        fastq_result = 0
    return fastq_result


verify_input_parameters( processor_number, parameters_file )
script_path = ""
reference_path = ""
with open( parameters_file ) as a:
    for line in a:
        if line.startswith( "#" ) is False:
            s = map( strip, line.split( ":" ) )
            if s[0] == "script_path":
                script_path = s[1]
                if os.path.exists( script_path ) is False:
                    print "The MetaShot path is inexistent"
                    sys.exit( )
            if s[0] == "reference_path":
                reference_path = s[1]
                if os.path.exists( reference_path ) is False:
                    print "The reference path is inexistent"
                    sys.exit( )

if script_path == "":
    print "The parameter file lacks of the script_path info"
    sys.exit( )
if reference_path == "":
    print "The parameter file lacks of the reference_path info"
    sys.exit( )

wd = os.getcwd( )
###################################
#    PRE-PROCESSING PROCEDURE     #
###################################
processor_number = 10
index = 1
data_processing_list = {}
# STDOUT and STDERR are redirected to the same file
std_out_file = open( "faqcs_stdout.out", "w" )
for data in multiple_input_data:
    data_processing_list.setdefault( index, [] )
    R1 = data[0]
    R2 = data[1]
    output_folder = os.path.join( wd, "trimmed_data_" + str( index ) )
    cmd = shlex.split( "FaQCs.pl -p %s %s -mode BWA_plus -q 25 -min_L 50 -n 2 -lc 0.70 -t %i  %s -d %s" % (
        R1, R2, processor_number, phix, output_folder) )
    print cmd
    p = subprocess.Popen( cmd, stdout=std_out_file, stderr=std_out_file )
    data_processing_list[index].append( p.pid )
    data_processing_list[index].append( output_folder )
    data_processing_list[index].append( R1 )
    data_processing_list[index].append( R2 )
    index += 1

completed = set( )
while len( completed ) != len( data_processing_list ):
    for index in data_processing_list.keys( ):
        if index not in completed:
            # print split
            pid = data_processing_list[index][0]
            exec_folder = data_processing_list[index][1]
            process_status = ""
            if psutil.pid_exists( pid ):
                process_status = psutil.Process( pid ).status( )
                # print pid, process_status
                # print len(completed), len(data_processing_list)
            if psutil.pid_exists( pid ) is False or process_status.lower( ) in ["finished", "zombie"]:
                trimmed_data_1 = os.path.join( exec_folder, "QC.1.trimmed.fastq" )
                trimmed_data_2 = os.path.join( exec_folder, "QC.2.trimmed.fastq" )
                result = verify_fastq( trimmed_data_1, trimmed_data_2 )
                if result == 0:
                    print exec_folder
                    R1 = data_processing_list[index][2]
                    R2 = data_processing_list[index][3]
                    del data_processing_list[index]
                    data_processing_list.setdefault( index, [] )
                    cmd = shlex.split( "FaQCs.pl -p %s %s -mode BWA_plus -q 25 -min_L 50 -n 2 -lc 0.70 -t %i -d %s" % (
                        R1, R2, processor_number, exec_folder) )
                    p = subprocess.Popen( cmd, stdout=std_out_file, stderr=std_out_file )
                    data_processing_list[index].append( p.pid )
                    data_processing_list[index].append( exec_folder )
                    data_processing_list[index].append( R1 )
                    data_processing_list[index].append( R2 )
                elif result >= 1:
                    completed.add( index )
std_out_file.close( )

###################################
#    FIND MICROBIAL CANDIDATES    #
###################################
for i in range( len( multiple_input_data ) ):
    i += 1
    folder = os.path.join( wd, "trimmed_data_%i" % i )
    os.chdir( folder )
    tmp = open( "read_list", "w" )
    tmp.write( "QC.1.trimmed.fastq\tQC.2.trimmed.fastq\n" )
    tmp.close( )
    cmd = shlex.split( "python %s -i read_list -r %s" % (os.path.join( script_path, "find_microbiome.py" ), reference_path) )
    p = subprocess.Popen( cmd )
    p.wait( )
    os.chdir( wd )

###################################
# CLEANED READ-LIST FILE CREATION #
###################################
# Annotazione dei file fastq contenenti le read denoised in un nuovo file read list
tmp = open( "read_list_cleaned", "w" )
for i in range( len( multiple_input_data ) ):
    i += 1
    for line in open( os.path.join( wd, "trimmed_data_%i" % i, "read_list" ) ):
        cleaned_list_file = map( strip, line.split( "\t" ) )
        tmp.write( "%s\t%s\n" % (os.path.join( wd, "trimmed_data_%i" % i, cleaned_list_file[0] ),
                               os.path.join( wd, "trimmed_data_%i" % i, cleaned_list_file[1] )) )
tmp.close( )

# mapping on the human genome and/or transcriptome
# to map all the human sequence against the human genome, first load the star reference genome
cmd = shlex.split( "STAR --genomeDir %s --genomeLoad LoadAndExit" % os.path.join(reference_path, "Homo_sapiens") )
p = subprocess.Popen( cmd )
p.wait( )
data_processing_list = {}
for i in range( len( multiple_input_data ) ):
    i += 1
    folder = os.path.join( wd, "trimmed_data_%i" % i )
    print folder
    data_processing_list.setdefault( i, [] )
    os.chdir( folder )
    mapper = shlex.split(
        "python %s -s human -i read_list -g -r %s" % (os.path.join( script_path, "host_mapper.py" ), reference_path) )
    p = subprocess.Popen( mapper )
    print p.pid
    data_processing_list[i].append( p.pid )
    data_processing_list[i].append( os.path.join( folder, "mapping_on_human", "human_dataAligned.out.sam" ) )
    os.chdir( wd )

completed = set( )
while len( completed ) != len( data_processing_list ):
    for index in data_processing_list.keys( ):
        if index not in completed:
            # print split
            pid = data_processing_list[index][0]
            sam_output = data_processing_list[index][1]
            process_status = ""
            if psutil.pid_exists( pid ):
                process_status = psutil.Process( pid ).status( )
                # print pid, process_status
                # print len(completed), len(data_processing_list)
            if psutil.pid_exists( pid ) is False or process_status.lower( ) in ["finished", "zombie"]:
                print sam_output
                result = control_mapping_procedure( sam_output )
                print result
                if result == 0:
                    folder = os.path.join( wd, "trimmed_data_%i" % index )
                    print folder
                    os.chdir( folder )
                    del data_processing_list[index]
                    data_processing_list.setdefault( index, [] )
                    mapper = shlex.split(
                        "python %s -s human -i read_list -g -r %s" % (os.path.join( script_path, "host_mapper.py" ), reference_path) )
                    p = subprocess.Popen( mapper )
                    data_processing_list[index].append( p.pid )
                    data_processing_list[index].append(
                        os.path.join( folder, "mapping_on_human", "human_dataAligned.out.sam" ) )
                    os.chdir( wd )
                elif result == 1:
                    completed.add( index )
                    print len( completed ), len( data_processing_list )

cmd = shlex.split( "STAR --genomeDir /home/bfosso/share/STAR/human --genomeLoad Remove" )
p = subprocess.Popen( cmd )
p.wait( )
human_folder = os.path.join( wd, "mapping_on_human" )
if os.path.exists( human_folder ) is False:
    os.mkdir( human_folder )
tmp = open( os.path.join( human_folder, "mapped_on_host.lst" ), "w" )
for i in completed:
    data = os.path.join( wd, "trimmed_data_%i" % i, "mapping_on_human", "mapped_on_host.lst" )
    with open( data ) as a:
        for line in a:
            tmp.write( line )
tmp.close( )

###################################
# MICROBIAL CANDIDATES READ LIST  #
###################################
# Annotazione dei file fastq contenenti le read denoised in un nuovo file read list
tmp = open( "candidate_microbial_list", "w" )
for i in range( len( multiple_input_data ) ):
    i += 1
    trimmed = os.path.join( wd, "trimmed_data_%i" % i )
    tmp.write( "%s\t%s\n" % (os.path.join( wd, "trimmed_data_%i" % i, "R1_micro_candidates.fastq" ),
                           os.path.join( wd, "trimmed_data_%i" % i, "R2_micro_candidates.fastq" )) )
tmp.close( )

###################################
#         DIVISION MAPPING        #
###################################
# procedure di mapping sulle divisioni
# Esecuzione sui batteri
bacterial_mapper = shlex.split(
    "python %s -i candidate_microbial_list -p 30 -s %s -r %s" % (os.path.join( script_path, "bacterial_division_anlyser.py" ), script_path, reference_path) )
p = subprocess.Popen( bacterial_mapper )
p.wait( )

multiple_division_process = {}
for division in ["virus", "fungi", "protist"]:
    multiple_division_process.setdefault( division, [] )
    division_mapper = shlex.split(
        "python %s -i candidate_microbial_list -d %s" % (os.path.join( script_path, "Division_analyser.py" ), division) )
    p = subprocess.Popen( division_mapper, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
    print p.pid
    multiple_division_process[division].append( p.pid )
    multiple_division_process[division].append( os.path.join( wd, division, "mapped_on_%s_total.txt" % division ) )

completed = set( )
while len( completed ) != len( multiple_division_process ):
    for division in multiple_division_process.keys( ):
        if division not in completed:
            # print split
            pid = multiple_division_process[division][0]
            mapped_reads_file = multiple_division_process[division][1]
            process_status = ""
            if psutil.pid_exists( pid ):
                process_status = psutil.Process( pid ).status( )
                # print pid, process_status
                # print len(completed), len(data_processing_list)
            if psutil.pid_exists( pid ) is False or process_status.lower( ) in ["finished", "zombie"]:
                if os.path.exists( mapped_reads_file ) and os.stat( mapped_reads_file )[6] != 0:
                    completed.add( division )
                else:
                    del multiple_division_process[division]
                    division_mapper = shlex.split(
                        "python %s -i candidate_microbial_list -d %s" % (os.path.join( script_path, "Division_analyser.py" ), division) )
                    p = subprocess.Popen( division_mapper, stdout=subprocess.PIPE, stderr=subprocess.STDOUT )
                    print p.pid
                    multiple_division_process.setdefault( division, [] )
                    multiple_division_process[division].append( p.pid )
                    multiple_division_process[division].append( os.path.join( wd, division, "mapped_on_%s_total.txt" % division ) )

###################################
#      PREPARE RESULT FILES       #
###################################
classifier = shlex.split(
    "python %s -s %s -r %s" % (os.path.join( script_path, "new_division_classifier.py" ), script_path , reference_path) )
p = subprocess.Popen( classifier )
p.wait( )

krona = shlex.split(
    "python %s " % (os.path.join( script_path, "krona_data_creator.py" )) )
p = subprocess.Popen( krona )
p.wait( )
