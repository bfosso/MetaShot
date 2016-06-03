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


def usage():
    print ('This script performs the taxonomic classification, based on the Division data\n'
           'Options'
           '\t-r\treference path [MANDATORY]\n'
           '\t-s\tscript path path [MANDATORY]\n'
           '\t-h    print this help.\n'
           'Usage:\n'
           '\tpython new_division_classifier.py -s host \n'
           '\t')


script_path = "/home/bfosso/share/MetaShot_reference/script/MetaShot_execution"

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


def define_correct_host_name(species):
    synonyms = {"human": ["Homo sapiens", "homo", "human"], "mouse": ["Mus musculus", "mouse", "mice", "mus"]}
    for KEY in synonyms.keys( ):
        if species.lower( ) in synonyms[KEY]:
            return KEY


host = define_correct_host_name( host )
if host is None:
    print "UNCORRECT HOST SPECIES"
    sys.exit( )

human_viruses = set( )
if host == "human":
    with open( os.path.join( reference_path, "Virus/human_viruses_selection/accepted_accession.lst" ) ) as a:
        for line in a:
            human_viruses.add( line.strip( ) )

division_list = [virus_set, "bacteria", "fungi", "protist"]
division_to_dump = dict( virus=os.path.join( reference_path, "Virus/tango_dmp_folder/VIRUS_DMP" ),
                         bacteria=os.path.join( reference_path, "Bacteria/tango_dmp_folder/BACTERIA_DMP" ),
                         fungi=os.path.join( reference_path, "Fungi/tango_dmp_folder/FUNGI_DMP" ),
                         protist=os.path.join( reference_path, "Protist/tango_dmp_folder/PROTISTA_DMP" ) )

division_to_visualization = dict( virus=os.path.join( reference_path, "Virus/tango_dmp_folder/visualization_virus.dmp" ),
                                  bacteria=os.path.join( reference_path, "Bacteria/tango_dmp_folder/visualization_bacterial.dmp" ),
                                  fungi=os.path.join( reference_path, "Fungi/tango_dmp_folder/visualization_fungi.dmp" ),
                                  protist=os.path.join( reference_path, "Protist/tango_dmp_folder/visualization_protista.dmp" ) )

herv2accepted_taxonomy = dict( virus=os.path.join( reference_path, "Virus/tango_dmp_folder/herv.lst" ) )
# creazione ambiguos list

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
host_folder = "mapping_on_" + host
if os.path.exists( host_folder ):
    print host
    for acc in open( os.path.join( host_folder, "mapped_on_host.lst" ) ):
        host_list.add( acc.strip( ) )
else:
    print "no", host_list, "data"
    sys.exit( )

tmp = open( "table_data", "w" )
print >> tmp, "Division\tMatch\tRefSeq_Virus\tBacteria\tFungi\tProtist\tHuman"
for key in division2set.keys( ):
    lista = [str( len( division2set[key] ) ), str( len( division2set[key].intersection( division2set[virus_set] ) ) ),
             str( len( division2set[key].intersection( division2set["bacteria"] ) ) ),
             str( len( division2set[key].intersection( division2set["fungi"] ) ) ),
             str( len( division2set[key].intersection( division2set["protist"] ) ) ),
             str( len( division2set[key].intersection( host_list ) ) )]
    print >> tmp, key + "\t" + "\t".join( lista )
tmp.close( )

if os.path.exists( "New_TANGO_perl_version" ) is False:
    shutil.copytree( os.path.join( reference_path, "New_TANGO_perl_version" ), "New_TANGO_perl_version" )

wd = os.getcwd( )
ambiguos_read = set( )
for division in division_list:
    investigable_division = list( division_list )
    investigable_division.remove( division )
    print investigable_division
    print division_list
    other_mapped_acc = set( ).union( host_list, division2set[investigable_division[0]], division2set[investigable_division[1]], division2set[investigable_division[2]] )
    if os.path.exists( division ):
        print division
        print "preparing data for taxonomic classification"
        tmp = open( os.path.join( wd, division, "only_%s_mapped_data.txt" % division ), "w" )
        NODESFILE = division_to_visualization[division]
        GB2tango_dmp = {}
        for line in open( NODESFILE ):
            fields = map( strip, line.split( "|" ) )
            s = fields[0].split( "#" )
            GB2tango_dmp[s[1]] = s[0]
        assigned = 0
        for line in open( os.path.join( wd, division, "mapped_on_%s_total.txt" % division ) ):
            s = map( strip, line.split( " " ) )
            acc = s[0]
            if acc in other_mapped_acc:
                ambiguos_read.add( acc )
            else:
                if len( human_viruses ) == 0:
                    stringa = [acc]
                    for i in s[1:]:
                        if GB2tango_dmp.has_key( i ):
                            stringa.append( GB2tango_dmp[i] )
                    if len( stringa ) > 1:
                        print >> tmp, " ".join( stringa )
                        assigned += 1
                    else:
                        ambiguos_read.add( acc )
                else:
                    stringa_human = [acc]
                    stringa_not_human = [acc]
                    for i in s[1:]:
                        if GB2tango_dmp.has_key( i ):
                            if i in human_viruses:
                                stringa_human.append( GB2tango_dmp[i] )
                            else:
                                stringa_not_human.append( GB2tango_dmp[i] )
                    if len( stringa_human ) > 1:
                        print >> tmp, " ".join( stringa_human )
                        assigned += 1
                    elif len( stringa_not_human ) > 1:
                        print >> tmp, " ".join( stringa_not_human )
                        assigned += 1
                    else:
                        ambiguos_read.add( acc )
        tmp.close( )
        if os.stat( os.path.join( wd, division, "only_%s_mapped_data.txt" % division ) )[6] != 0:
            os.chdir( "New_TANGO_perl_version" )
            dump = division_to_dump[division]
            match_file = os.path.join( wd, division, "only_%s_mapped_data.txt" % division )
            output_file = os.path.join( wd, division, "%s_classification_data.txt" % division )
            cmd = shlex.split(
                "perl tango.pl --taxonomy %s --matches %s --output %s" % (dump, match_file, output_file) )
            # p = subprocess.Popen( cmd )
            # p.wait( )
            os.chdir( wd )
        else:
            print "no classificable sequences for " + division
    else:
        print "no data for " + division

tmp = open( os.path.join( host_folder, "only_%s.lst" % host ), "w" )
for acc in open( os.path.join( host_folder, "mapped_on_host.lst" ) ):
    if acc.strip( ) not in ambiguos_read:
        print >> tmp, acc.strip( )
tmp.close( )

if host == "human":
    accepted = division2set[virus_set].intersection( host_list )
    print len( accepted )
    exclusion = set( ).union( division2set["bacteria"], division2set["fungi"], division2set["protist"] )
    print len( exclusion )
    prob_herv = accepted.difference( exclusion )
    print len( prob_herv )
    if len( prob_herv ) > 0:
        print "PROB HERV", len( prob_herv )
        GB2tango_dmp = {}
        for line in open( division_to_visualization[virus_set] ):
            line = line.strip( )
            fields = map( strip, line.split( "|" ) )
            s = fields[0].split( "#" )
            GB2tango_dmp[s[1]] = s[0]
        herv_accepted_taxonomy = set( )
        for line in open( herv2accepted_taxonomy[virus_set] ):
            s = map( strip, line.split( "\t" ) )
            herv_accepted_taxonomy.add( s[0] )
        folder = os.path.join( wd, "human_retroviruses" )
        os.mkdir( folder )
        tmp = open( os.path.join( folder, "only_human_HERV_mapped_data.txt" ), "w" )
        for line in open( os.path.join( virus_set, "mapped_on_%s_total.txt" % virus_set ) ):
            s = map( strip, line.split( " " ) )
            if s[0] in prob_herv:
                stringa = [s[0]]
                for i in s[1:]:
                    if GB2tango_dmp.has_key( i ):
                        stringa.append( GB2tango_dmp[i] )
                if len( stringa ) > 1:
                    print >> tmp, " ".join( stringa )
        tmp.close( )
        if os.stat( os.path.join( folder, "only_human_HERV_mapped_data.txt" ) )[6] != 0:
            os.chdir( "New_TANGO_perl_version" )
            dump = division_to_dump[virus_set]
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
            ambiguos_read = ambiguos_read.difference( found_herv )
            with open( os.path.join( wd, virus_set, "%s_classification_data.txt_0.5" % virus_set ), "a" ) as viral_ass:
                for line in open( final_output ):
                    viral_ass.write( line )
    os.chdir( wd )
    lista = list( ambiguos_read )
    tmp = open( "ambiguos_pe_read.lst", "w" )
    tmp.write( "\n".join( lista ) )
    tmp.close( )
else:
    lista = list( ambiguos_read )
    tmp = open( "ambiguos_pe_read.lst", "w" )
    tmp.write( "\n".join( lista ) )
    tmp.close( )

for division in division_list:
    output_file = os.path.join( wd, division, "%s_classification_data.txt" % division )
    if os.path.exists( output_file + "_0.5" ) and os.stat( output_file + "_0.5" )[6] != 0:
        os.chdir( os.path.join( wd, division ) )
        NODESFILE = division_to_visualization[division]
        cmd = shlex.split(
            "python %s -n %s -d %s -i %s" % (os.path.join( script_path, "tree_builder_for_perl_tango.py" ), NODESFILE, division, output_file) )
        p = subprocess.Popen( cmd )
        p.wait( )
        os.chdir( wd )
        tree_name = os.path.join( wd, division, division + "_tree.nwk" )
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
        os.chdir( wd )
