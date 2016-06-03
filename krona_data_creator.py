__author__ = 'Bruno Fosso'
__version__ = "1.0"

import getopt
import os
import shlex
import subprocess
from string import strip


def usage():
    print ('This script performs the taxonomic classification, based on the Division data\n'
           'Options'

           '\t-v    Virus representative set. Available sets are virus and refseq_virus [DEFAULT virus]\n'
           '\t-h    print this help.\n'
           'Usage:\n'
           '\tpython new_division_classifier.py -s host \n'
           '\t')


reference_path = ""
try:
    import sys

    opts, args = getopt.getopt( sys.argv[1:], "r:v:h" )
except getopt.GetoptError, err:
    print str( err )
    usage( )
    sys.exit( )
for o, a in opts:
    if o == "-h":
        usage( )
        sys.exit( )
    elif o == "-r":
        reference_path = a
    elif o == "-v":
        virus_set = a
        if virus_set not in ["virus", "refseq_virus"]:
            print usage( )
            exit( )
    else:
        usage( )
        exit( )

if reference_path == "":
    print "The parameter file lacks of the reference_path info"
    sys.exit( )

host = "human"
virus_set = "virus"


def define_correct_host_name(species):
    synonyms = {"human": ["Homo sapiens", "homo", "human"], "mouse": ["Mus musculus", "mouse", "mice", "mus"]}
    for key in synonyms.keys( ):
        if species.lower( ) in synonyms[key]:
            return key


host = define_correct_host_name( host )
if host is None:
    print "UNCORRECT HOST SPECIES"
    sys.exit( )

if host == "human":
    human_viruses = set( )
    host_taxid = "9606"
    for line in open( "/home/bfosso/share/MetaShot_reference/Virus/human_viruses_selection/accepted_accession.lst" ):
        human_viruses.add( line.strip( ) )
else:
    host_taxid = "10090"
division_list = [virus_set, "bacteria", "fungi", "protist"]

division_to_visualization = dict( virus=os.path.join( reference_path, "Virus/tango_dmp_folder/visualization_virus.dmp" ),
                                  bacteria=os.path.join( reference_path, "Bacteria/tango_dmp_folder/visualization_bacterial.dmp" ),
                                  fungi=os.path.join( reference_path, "Fungi/tango_dmp_folder/visualization_fungi.dmp" ),
                                  protist=os.path.join( reference_path, "Protist/tango_dmp_folder/visualization_protista.dmp" ) )

krona_folder = os.path.join( os.getcwd( ), "krona_data" )
if os.path.exists( krona_folder ) is not True:
    os.mkdir( krona_folder )

total = open( os.path.join( krona_folder, "full_data_krona_data.tsv" ), "w" )
for division in division_list:
    print division
    NODESFILE = division_to_visualization[division]
    node2name = {}
    node2order = {}
    for line in open( NODESFILE ):
        line = line.strip( )
        fields = map( strip, line.split( "|" ) )
        s = fields[0].split( "#" )
        nodeid = s[0]
        node2name[s[0]] = s[1]
        node2order[s[0]] = s[2]
    if os.path.exists( os.path.join( division, division + "_classification_data.txt_0.5" ) ):
        tmp = open( os.path.join( krona_folder, division + "_krona_data.tsv" ), "w" )
        with open( os.path.join( division, division + "_classification_data.txt_0.5" ) ) as a:
            for line in a:
                s = map( strip, line.split( "\t" ) )
                read, path = s[0], s[2].split( ";" )
                order = node2order[path[0]]
                if order == "GB acc":
                    node = path[1]
                else:
                    node = path[0]
                # print node2order[node]
                print >> tmp, read + "\t" + node
                print >> total, read + "\t" + node
        tmp.close( )
        print "krona generation"
        cmd = shlex.split( "ktImportTaxonomy -o %s -k %s -tax %s" % (
            division + "_krona.html", os.path.join( krona_folder, division + "_krona_data.tsv" ), os.path.join( reference_path, "krona_tax" )) )
        p = subprocess.Popen( cmd )
        p.wait( )
    else:
        print "no classification data for %s" % division

host_data = os.path.join( "mapping_on_" + host, "only_human.lst" )
with open( host_data ) as a:
    for line in a:
        print >> total, line.strip( ) + "\t" + host_taxid
total.close( )
cmd = shlex.split( "ktImportTaxonomy -o %s -k %s -tax %s" % (
    "full_assignment_krona.html", os.path.join( krona_folder, "full_data_krona_data.tsv" ), os.path.join( reference_path, "krona_tax" )) )
p = subprocess.Popen( cmd )
p.wait( )
