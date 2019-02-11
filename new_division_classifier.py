__author__ = 'Bruno Fosso'
__version__ = "1.0"

import os
import subprocess
import shlex
from string import strip
from ete2 import Tree
import shutil
import getopt
import sys
import psutil
import numpy
from collections import OrderedDict
from operator import itemgetter


def usage():
    print ('This script performs the taxonomic classification, based on the Division data\n'
           'Options'
           '\t-r\treference path [MANDATORY]\n'
           '\t-s\tscript path path [MANDATORY]\n'
           '\t-h    print this help.\n'
           'Usage:\n'
           '\tpython new_division_classifier.py -s /path/to/script -r /path/to/reference \n'
           '\t')


script_path = ""
reference_path = ""
try:
    import sys

    opts, args = getopt.getopt( sys.argv[1:], "s:v:r:h" )
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
    elif o == "-s":
        script_path = a
    elif o == "-v":
        virus_set = a
        if virus_set not in ["virus", "refseq_virus"]:
            print usage( )
            exit( )
    else:
        usage( )
        exit( )

if script_path == "":
    print "The parameter file lacks of the script_path info"
    sys.exit( )
if reference_path == "":
    print "The parameter file lacks of the reference_path info"
    sys.exit( )

virus_set = "virus"
host = "human"


def define_correct_host_name(host_species):
    """
    Control the corrct name for host species
    :param host_species:
    :return: A key value
    """
    synonyms = {"human": ["Homo sapiens", "homo", "human"]}
    for KEY in synonyms.keys( ):
        if host_species.lower( ) in synonyms[KEY]:
            return KEY


host = define_correct_host_name( host )
if host is None:
    print "UNCORRECT HOST SPECIES"
    sys.exit( )

human_viruses = set( )
with open( os.path.join( reference_path, "Virus/human_viruses_selection/accepted_accession.lst" ) ) as a:
    for line in a:
        human_viruses.add( line.strip( ) )

division_list = [virus_set, "bacteria", "fungi", "protist"]

acc2node = {}
acc2division = {}
with open( os.path.join( reference_path, "Metashot_reference_taxonomy/ACC_num_to_division_info.tsv" ) ) as a:
    for line in a:
        field = map( strip, line.split( "\t" ) )
        acc, node, ref_division = field[0], field[1], field[2]
        acc2node[acc] = node
        acc2division[acc] = ref_division

dump = os.path.join( reference_path, "Metashot_reference_taxonomy/cleaned_ncbi_ref.bin" )

herv2accepted_taxonomy = dict( virus=os.path.join( reference_path, "Virus/tango_dmp_folder/herv.lst" ) )
# creazione ambiguos list


node2parent = {}
node2name = {}
node2order = {}
with open( os.path.join( reference_path, "Metashot_reference_taxonomy/nodes.dmp" ) ) as nodes:
    for line in nodes:
        fields = map( strip, line.split( "|" ) )
        node, parent, order = fields[0], fields[1], fields[2]
        node2parent[node] = parent
        node2order[node] = order

with open( os.path.join( reference_path, "Metashot_reference_taxonomy/names.dmp" ) ) as names:
    for line in names:
        fields = map( strip, line.split( "|" ) )
        if "scientific name" in fields:
            node, name = fields[0], fields[1]
            node2name[node] = name

virus_list = set( )
bacteria_list = set( )
fungi_list = set( )
protist_list = set( )
division2set = {virus_set: virus_list, "bacteria": bacteria_list, "fungi": fungi_list, "protist": protist_list}

for division in division_list:
    if os.path.exists( division ):
        print division
        with open( os.path.join( division, "mapped_on_%s_total.txt" % division ) ) as c:
            for line in c:
                acc = line.split( " " )[0]
                division2set[division].add( acc )
    else:
        print "no", division, "data"

host_list = set( )
host_folder = "mapping_on_%s" % host
if os.path.exists( host_folder ):
    print host
    for acc in open( os.path.join( host_folder, "mapped_on_host.lst" ) ):
        host_list.add( acc.strip( ) )
else:
    print "no", host_list, "data"
    sys.exit( )

tmp = open( "table_data", "w" )
print >> tmp, "Division\tMatch\tRefSeq_Virus\tBacteria\tFungi\tProtist\tHuman"
all_mapped_reads = host_list
for key in division2set.keys( ):
    all_mapped_reads = all_mapped_reads.union( division2set[key] )
    lista = [str( len( division2set[key] ) ), str( len( division2set[key].intersection( division2set[virus_set] ) ) ),
             str( len( division2set[key].intersection( division2set["bacteria"] ) ) ),
             str( len( division2set[key].intersection( division2set["fungi"] ) ) ),
             str( len( division2set[key].intersection( division2set["protist"] ) ) ),
             str( len( division2set[key].intersection( host_list ) ) )]
    print >> tmp, key + "\t" + "\t".join( lista )
tmp.close( )

ambiguos_read = {}
for reads in all_mapped_reads:
    found = []
    if reads in host_list:
        found.append( "human" )
    if reads in division2set[virus_set]:
        found.append( "virus" )
    if reads in division2set["bacteria"]:
        found.append( "bacteria" )
    if reads in division2set["fungi"]:
        found.append( "fungi" )
    if reads in division2set["protist"]:
        found.append( "protist" )
    if len( found ) > 1:
        ambiguos_read[reads] = found

tmp = open( os.path.join( host_folder, "only_%s.lst" % host ), "w" )
with open( os.path.join( host_folder, "mapped_on_host.lst" ) ) as human_mapped_data:
    for line in human_mapped_data:
        if ambiguos_read.has_key( line.strip( ) ) is False:
            print >> tmp, line.strip( )
tmp.close( )

if os.path.exists( "New_TANGO_perl_version" ) is False:
    shutil.copytree( os.path.join( reference_path, "New_TANGO_perl_version" ), "New_TANGO_perl_version" )

working_dir = os.getcwd( )
for division in division_list:
    if os.path.exists( division ):
        print division
        print "preparing data for taxonomic classification"
        tmp = open( os.path.join( working_dir, division, "only_%s_mapped_data.txt" % division ), "w" )
        assigned = 0
        with open( os.path.join( working_dir, division, "mapped_on_%s_total.txt" % division ) ) as md:
            for line in md:
                mapping_result = map( strip, line.split( " " ) )
                acc = mapping_result[0]
                if ambiguos_read.has_key( acc ) is False:
                    printing_string = [acc]
                    stringa_human = set( )
                    stringa_not_human = set( )
                    for i in mapping_result[1:]:
                        if i in human_viruses:
                            stringa_human.add( acc2node[i] )
                        else:
                            stringa_not_human.add( acc2node[i] )
                    if len( stringa_human ) > 0:
                        print >> tmp, " ".join( printing_string + list( stringa_human ) )
                        assigned += 1
                    elif len( stringa_not_human ) > 0:
                        print >> tmp, " ".join( printing_string + list( stringa_not_human ) )
                        assigned += 1
        tmp.close( )
        print division,assigned
        if os.stat( os.path.join( working_dir, division, "only_%s_mapped_data.txt" % division ) )[6] != 0:
            os.chdir( "New_TANGO_perl_version" )

            match_file = os.path.join( working_dir, division, "only_%s_mapped_data.txt" % division )
            output_file = os.path.join( working_dir, division, "%s_classification_data.txt" % division )
            cmd = shlex.split(
                "perl tango.pl --taxonomy %s --matches %s --output %s" % (dump, match_file, output_file) )
            p = subprocess.Popen( cmd )
            p.wait( )
            os.chdir( working_dir )
        else:
            print "no classificable sequences for %s" % division
    else:
        print "no data for %s" % division


if host == "human":
    accepted = division2set[virus_set].intersection( host_list )
    # print len( accepted )
    exclusion = set( ).union( division2set["bacteria"], division2set["fungi"], division2set["protist"] )
    # print len( exclusion )
    prob_herv = accepted.difference( exclusion )
    # print len( prob_herv )
    if len( prob_herv ) > 0:
        print "PROB HERV", len( prob_herv )
        herv_accepted_taxonomy = set( )
        for line in open( herv2accepted_taxonomy[virus_set] ):
            s = map( strip, line.split( "\t" ) )
            herv_accepted_taxonomy.add( s[0] )
        folder = os.path.join( working_dir, "human_retroviruses" )
        if os.path.exists( folder ) is False:
            os.mkdir( folder )
        tmp = open( os.path.join( folder, "only_human_HERV_mapped_data.txt" ), "w" )
        for line in open( os.path.join( virus_set, "mapped_on_%s_total.txt" % virus_set ) ):
            s = map( strip, line.split( " " ) )
            if s[0] in prob_herv:
                stringa = [s[0]]
                for i in s[1:]:
                    if acc2node.has_key( i ):
                        stringa.append( acc2node[i] )
                if len( stringa ) > 1:
                    print >> tmp, " ".join( stringa )
        tmp.close( )
        if os.stat( os.path.join( folder, "only_human_HERV_mapped_data.txt" ) )[6] != 0:
            os.chdir( "New_TANGO_perl_version" )
            match_file = os.path.join( folder, "only_human_HERV_mapped_data.txt" )
            output_file = os.path.join( folder, "human_HERV_classification_data.txt" )
            cmd = shlex.split(
                "perl tango.pl --taxonomy %s --matches %s --output %s" % (dump, match_file, output_file) )
            p = subprocess.Popen( cmd )
            p.wait( )
            found_herv = set( )
            final_output = os.path.join( folder, "human_HERV.txt" )
            if os.path.exists( output_file + "_0.5" ) and os.stat( output_file + "_0.5" )[6] != 0:
                tmp = open( final_output, "w" )
                for line in open( output_file + "_0.5" ):
                    s = map( strip, line.split( "\t" ) )
                    path = set( s[2].split( ";" ) )
                    if len( path.intersection( herv_accepted_taxonomy ) ) >= 1:
                        tmp.write( line )
                        found_herv.add( s[0] )
                tmp.close( )
            print "HERV", len( found_herv )
            for read in found_herv:
                ambiguos_read.pop( read, None )
            with open( os.path.join( working_dir, virus_set, "%s_classification_data.txt_0.5" % virus_set ), "a" ) as viral_ass:
                with open( final_output ) as herv_result:
                    for line in herv_result:
                        viral_ass.write( line )
    os.chdir( working_dir )
    tmp = open( "ambiguos_pe_read.lst", "w" )
    for read in ambiguos_read:
        tmp.write( "%s\n" % "\t".join( [read] + ambiguos_read[read] ) )
    tmp.close( )
else:
    tmp = open( "ambiguos_pe_read.lst", "w" )
    for read in ambiguos_read:
        tmp.write( "%s\n" % "\t".join( [read] + ambiguos_read[read] ) )
    tmp.close( )

####################
#  TAXON REFINING  #
####################

species2reads = {}
read2macth_list = {}
tot_species_assingment = 0
read2assigned_node = {}
for division in division_list:
    if os.path.exists( os.path.join( working_dir, division, "%s_classification_data.txt_0.5" % division ) ):
        with open( os.path.join( working_dir, division, "%s_classification_data.txt_0.5" % division ) ) as tango_assignment_path:
            for line in tango_assignment_path:
                s = map( strip, line.split( "\t" ) )
                read = s[0]
                path = s[2].split( ";" )
                for node in path:
                    if node2order.has_key( node ) and node2order[node] == "species":
                        species2reads.setdefault( node, set( ) )
                        species2reads[node].add( read )
                        tot_species_assingment += 1
                        read2assigned_node[read] = node
        with open( os.path.join( working_dir, division, "only_%s_mapped_data.txt" % division ) ) as mapping_data:
            for line in mapping_data:
                s = map( strip, line.split( " " ) )
                read = s[0]
                read2macth_list[read] = s[1:]

phix_path = []
node = "374840"
parent = node2parent[node]
while node != parent:
    phix_path.append( node )
    node = parent
    parent = node2parent[node]

for name in os.listdir( working_dir ):
    if name.startswith( "phix_removal_" ):
        if os.path.exists( os.path.join( working_dir, name, "Phix_sequences.lst" ) ):
            with open( os.path.join( working_dir, name, "Phix_sequences.lst" ) ) as phix_assignment:
                for line in phix_assignment:
                    species2reads.setdefault( "374840", set( ) )
                    species2reads["374840"].add( line.strip( ) )
                    tot_species_assingment += 1
                    read2macth_list[line.strip( )] = ["374840"]
                    read2assigned_node[line.strip( )] = "374840"

herv_read_list = set( )
print tot_species_assingment
if os.path.exists( os.path.join( working_dir, "human_retroviruses", "human_HERV.txt" ) ):
    with open( os.path.join( working_dir, "human_retroviruses", "human_HERV.txt" ) ) as herv_match:
        for line in herv_match:
            s = map( strip, line.split( "\t" ) )
            herv_read_list.add( s[0] )

species2remap = set( )
species2abundances = {}
for species in species2reads.keys( ):
    if (len( species2reads[species] ) / float(tot_species_assingment)) <= 0.00001:
        print node2name[species], len( species2reads[species] ) / tot_species_assingment
        species2remap.add( species )
        species2abundances[species] = len( species2reads[species] ) / tot_species_assingment
    else:
        species2abundances[species] = len( species2reads[species] ) / tot_species_assingment

for low_ab_species in species2remap:
    to_remove_reads = set( )
    for read in species2reads[low_ab_species]:
        if read not in herv_read_list:
            possible_species = set( )
            for node in read2macth_list[read]:
                parent = node2parent[node]
                while node != parent:
                    if node2order[node] == "species":
                        possible_species.add( node )
                        parent = node
                    else:
                        node = parent
                        parent = node2parent[node]
            abundances = {}
            for node in possible_species:
                if node != low_ab_species and species2abundances.has_key( node ) and node not in species2remap:
                    abundances[node] = species2abundances[node]
            if len( abundances ) != 0:
                abundances_ordered = OrderedDict( sorted( abundances.items( ), key=itemgetter( 1 ), reverse=True ) )
                read2assigned_node[read] = abundances.keys( )[0]
            else:
                read2assigned_node[read] = node2parent[read2assigned_node[read]]

for division in division_list:
    if os.path.exists( os.path.join( working_dir, division, "%s_classification_data.txt_0.5" % division ) ):
        new_taxonomic_assignment = open( os.path.join( working_dir, division, "%s_refined_classification_data.txt" % division ), "w" )
        with open( os.path.join( working_dir, division, "%s_classification_data.txt_0.5" % division ) ) as tango_assignment_path:
            for line in tango_assignment_path:
                s = map( strip, line.split( "\t" ) )
                if read2assigned_node.has_key(s[0]) and s[0] not in herv_read_list:
                    new_taxonomic_assignment.write( "%s\t%s\n" % (s[0], read2assigned_node[s[0]]) )
                elif s[0] in herv_read_list:
                    new_taxonomic_assignment.write( "%s\t%s\n" % (s[0], s[2].split(";")[0]) )
                else:
                    new_taxonomic_assignment.write( "%s\t%s\n" % (s[0], s[2].split( ";" )[0]) )

        new_taxonomic_assignment.close( )

for division in division_list:
    output_file = os.path.join( working_dir, division, "%s_refined_classification_data.txt" % division )
    if os.path.exists( output_file ) and os.stat( output_file )[6] != 0:
        os.chdir( os.path.join( working_dir, division ) )
        cmd = shlex.split(
            "python %s -d %s -i %s -n %s" % (os.path.join( script_path, "tree_builder_for_perl_tango.py" ), division, output_file, reference_path) )
        p = subprocess.Popen( cmd )
        p.wait( )
        os.chdir( working_dir )
        tree_name = os.path.join( working_dir, division, division + "_tree.nwk" )
        csv = open( division + "_CSV_result.html", "w" )
        csv.write( "Taxon Name\tRank\tTaxonomy ID\tDirectly Assigned Sequences\tDescendants Assignments\n" )
        if os.path.exists( tree_name ) and os.stat( tree_name )[6] != 0:
            t = Tree( tree_name )
            tree_testo = ""
            assigned = 0
            for node in t.iter_search_nodes( ):
                # print node.name
                if node.name == "NoName":
                    pass
                else:
                    assigned += int( node.assigned )
                    tree_testo += "['%s','%s','%s', %s, %s]," % (node.name.replace( "'", "" ), node.Order, node.taxid, node.assigned, node.summarized)
                    csv.write( "%s\t%s\t%s\t%s\t%s\n" % (node.name.replace( "'", "" ), node.Order, node.taxid, node.assigned, node.summarized) )
            csv.close( )
            tree_testo = tree_testo.rstrip( "," )
            if tree_testo != "":
                print "prepare graphical representation"
                stringa = """
                <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
                <html xmlns="http://www.w3.org/1999/xhtml">
                    <head>
                        <meta http-equiv="content-type" content="text/html; charset=utf-8"/>
                        <title>
                            %s
                        </title>
                        <script type="text/javascript" src="http://www.google.com/jsapi"></script>
                        <script type="text/javascript">
                            google.load('visualization', '1.1', {packages: ['controls']});
                            </script>
                        <script type="text/javascript">
                            function drawVisualization() {
                                // Prepare the data
                                var data = google.visualization.arrayToDataTable([['Taxon Name','Rank','Taxonomy ID','Directly Assigned Sequences','Descendants Assignments'],%s
                                ]);

                                // Define category pickers for 'Country', 'Region/State' and 'City'
                                var countryPicker = new google.visualization.ControlWrapper({
                                                                                            'controlType': 'CategoryFilter',
                                                                                            'containerId': 'control1',
                                                                                            'options': {
                                                                                            'filterColumnLabel': 'Rank',
                                                                                            'ui': {
                                                                                            'labelStacking': 'vertical',
                                                                                            'allowTyping': true,
                                                                                            'allowMultiple': true,},
                                                                                            },
                                                                                            'state': {'selectedValues':['Family']}
                                                                                            });

                                var regionPicker = new google.visualization.ControlWrapper({
                                                                                           'controlType': 'CategoryFilter',
                                                                                           'containerId': 'control2',
                                                                                           'options': {
                                                                                           'filterColumnLabel': 'Taxon Name',
                                                                                           'ui': {
                                                                                           'labelStacking': 'vertical',
                                                                                           'allowTyping': true,
                                                                                           'allowMultiple': false
                                                                                           }
                                                                                           }
                                                                                           });

                                var table = new google.visualization.ChartWrapper({
                                                                                  'chartType': 'Table',
                                                                                  'containerId': 'chart1',
                                                                                  'options': {'width': '1000px','height':'1000px'}
                                                                                  });
                                new google.visualization.Dashboard(document.getElementById('dashboard')).
                                bind(countryPicker,regionPicker).
                                bind(regionPicker,table).
                                // Draw the dashboard
                                draw(data);
                            }


                            google.setOnLoadCallback(drawVisualization);
                            </script>
                            </head>
                    <body style="font-family: Arial;border: 0 none;">
                    <h1 align="center"><em>Taxonomic Assignment Table for %s</em></h1>
                    <h2>Assigned sequences: %i</h2>
                        <div id="dashboard">
                            <table>
                                <tr style='vertical-align: top'>
                                    <td style='width: 300px; font-size: 0.9em;'>
                                        <div id="control1"></div>
                                        <div id="control2"></div>
                                    </td>
                                    <td style='width: 600px'>
                                        <div style="float: left;" id="chart1"></div>
                                </tr>
                            </table>
                        </div>
                    </body>
                </html>
                """ % (division.upper( ), tree_testo, division.upper( ), assigned)
                tmp = open( division + "_html_result.html", "w" )
                tmp.write( stringa )
                tmp.close( )
            else:
                print "no textual data tree for " + division
        else:
            print "no available tree for " + division
        os.chdir( working_dir )
