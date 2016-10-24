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
           '\t-r    reference path [MANDATORY].\n'
           '\t-h    print this help.\n'
           'Usage:\n'
           '\tpython new_division_classifier.py  \n'
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
    else:
        usage( )
        exit( )

if reference_path == "":
    print "The parameter file lacks of the reference_path info"
    sys.exit( )

host = "human"
virus_set = "virus"

def define_correct_host_name(species):
    synonyms = {"human": ["Homo sapiens", "homo", "human"]}
    for key in synonyms.keys( ):
        if species.lower( ) in synonyms[key]:
            return key


host = define_correct_host_name( host )
if host is None:
    print "UNCORRECT HOST SPECIES"
    sys.exit( )

##################################
# COSTRUZIONE DEI DIZIONARI	     #
##################################
node2parent = {}
node2name = {}
node2order = {}
with open(os.path.join(reference_path,"Metashot_reference_taxonomy/nodes.dmp")) as nodes:
    for line in nodes:
        fields = map(strip, line.split("|"))
        node, parent, order =  fields[0], fields[1], fields[2]
        node2parent[node] = parent
        node2order[node] = order

with open( os.path.join( reference_path, "Metashot_reference_taxonomy/names.dmp" ) ) as names:
    for line in names:
        fields = map( strip, line.split( "|" ) )
        if "scientific name" in fields:
            node, name = fields[0] ,fields[1]
            node2name[node] = name

if host == "human":
    human_viruses = set( )
    host_taxid = "9606"
    for line in open( "/home/bfosso/share/MetaShot_reference/Virus/human_viruses_selection/accepted_accession.lst" ):
        human_viruses.add( line.strip( ) )
else:
    host_taxid = "10090"
division_list = [virus_set, "bacteria", "fungi", "protist"]


krona_folder = os.path.join( os.getcwd( ), "krona_data" )
if os.path.exists( krona_folder ) is not True:
    os.mkdir( krona_folder )

total = open( os.path.join( krona_folder, "full_data_krona_data.tsv" ), "w" )
for division in division_list:
    print division
    if os.path.exists( os.path.join( division, "%s_refined_classification_data.txt" % division) ):
        tmp = open( os.path.join( krona_folder, "%s_krona_data.tsv" % division ), "w" )
        with open( os.path.join( division, "%s_refined_classification_data.txt" % division ) ) as a:
            for line in a:
                s = map( strip, line.split( "\t" ) )
                read, node = s[0], s[1]
                tmp.write("%s\t%s\n" % (read, node))
                total.write("%s\t%s\n" % (read, node))
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
