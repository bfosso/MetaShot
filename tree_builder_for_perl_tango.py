__author__ = 'Bruno Fosso'
import getopt
import os
import sys
from string import strip
from ete2 import Tree, TreeNode

def usage():
    print ('This script converts the taxonomic assignments made by TANGO in a taxonomical tree:\n'
           'Options:\n'
           '\t-n    reference path \n'
           '\t-d    division'
           '\t-i    taxonomic assignment'
           '\t-h    print this help.\n'
           'Usage:\n'
           '\tpython tango_ass_2_taxa_freq.py -n nodes_file -d division -i taxonomic assignment file\n'
           '\t'
    )


reference_path = ""
basename = ""
tango_output = ""
try:
    opts, args = getopt.getopt(sys.argv[1:], "hd:n:i:")
except getopt.GetoptError, err:
    print str(err)
    usage()
    sys.exit()
for o, a in opts:
    if o == "-h":
        usage()
        sys.exit()
    elif o == "-d":
        basename = a
    elif o == "-i":
        tango_output = a
    elif o == "-n":
        reference_path = a
    else:
        assert False, "Unhandled option."

if basename != "":
    print basename
else:
    print "no division"
    exit()


def controllo_tango_exec(l):
    """


    """
    if os.path.exists(l):
        pass
    else:
        print "no correct tango output"
        print "please read the tango_error.log file"
        sys.exit()


controllo_tango_exec(tango_output)

##################################
# COSTRUZIONE DEI DIZIONARI	     #
##################################
id2node = {}
node2parent = {}
node2name = {}
node2order = {}
all_ids = set([])
all_nodes = []
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

# DEPRECATED
#####################################
# ANALISI DEI DATI DAL TANGO OUTPUT #
#####################################
#
# out_acc = set()
# seq2taxa = {}
# result = open(tango_output)
# for line in result.readlines():
#     line = line.strip()
#     s = line.split("\t")
#     out_acc.add(s[0])
#     path = s[2].split(";")
#     if node2order[path[0]] == "GB acc":
#         seq2taxa[s[0]] = node2parent[path[0]]
#     else:
#         seq2taxa[s[0]] = path[0]
# result.close()
#####################################

#########################################
# ANALISI DEI DATI DAL TAXON REFINEMENT #
#########################################

out_acc = set()
seq2taxa = {}
with open(tango_output) as result:
    for line in result:
        ass_data = map(strip, line.split("\t"))
        seq, assignment_taxa = ass_data[0], ass_data[1]
        seq2taxa[seq] = assignment_taxa
        out_acc.add(seq)

taxa = set()
assigned = {}
for key in seq2taxa.keys():
    taxa.add(seq2taxa[key])
    s = key.split(";")
    if len(s) >= 2:
        value = int(s[1].lstrip("size="))
    else:
        value = 1
    if assigned.has_key(seq2taxa[key]):
        assigned[seq2taxa[key]] += value
    else:
        assigned[seq2taxa[key]] = value

##################################
# COSTRUZIONE ALBERO			 #
##################################
ass_node2parent = {"1": "1"}
for nodeid in taxa:
    parentid = node2parent[nodeid]
    while nodeid != parentid:  #costruiamo un nuovo dizionario per i soli taxa che abbiamo identificato nel campione
        ass_node2parent[nodeid] = parentid
        nodeid = parentid
        parentid = node2parent[nodeid]

node2parentid = {}
for nodeid in ass_node2parent.keys():
    parentid = ass_node2parent[nodeid]
    # Stores node connections
    all_ids.update([nodeid, parentid])
    # Creates a new TreeNode instance for each new node in file
    n = TreeNode()
    # Sets some TreeNode attributes
    n.add_feature("name", node2name[nodeid])
    n.add_feature("taxid", nodeid)
    n.add_feature("Order", node2order[nodeid])

    # updates node list and connections
    node2parentid[n] = parentid
    id2node[nodeid] = n
print len(id2node)
# Reconstruct tree topology from previously stored tree connections
print 'Reconstructing tree topology...'
for node in id2node.itervalues():
    parentid = node2parentid[node]
    parent = id2node[parentid]
    # node with taxid=1 is the root of the tree
    if node.taxid == "1":
        t = node
    else:
        parent.add_child(node)

freq = {}
for node in t.iter_search_nodes():
    if assigned.has_key(node.taxid):
        val = int(assigned[node.taxid])
        node.add_feature("assigned", assigned[node.taxid])
    else:
        val = 0
        node.add_feature("assigned", "0")
    for child in node.iter_descendants():
        if assigned.has_key(child.taxid):
            val = val + assigned[child.taxid]
    node.add_feature("summarized", str(val))

for node in t.iter_search_nodes(name="NoName"):
    if node.is_root():
        node.add_feature("assigned", "0")
        node.add_feature("Order", "root")
        count = 0
        for nodo in node.iter_descendants():
            count += int(nodo.assigned)
        node.add_feature("summarized", str(count))

open(basename + "_tree.nwk", "w").write(t.write(features=["name", "taxid", "assigned", "summarized", "Order"]))

t = Tree(basename + "_tree.nwk")


