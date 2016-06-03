__author__ = 'Bruno Fosso'
import getopt
import os
import sys
from string import strip
from ete2 import Tree, TreeNode

def usage():
    print ('This script converts the taxonomic assignments made by TANGO in a taxonomical tree:\n'
           'Options:\n'
           '\t-n    nodes file of the taxonomic reference tree. If not indicated it uses the bacterial taxonomy \n'
           '\t-d    division'
           '\t-i    tango output'
           '\t-h    print this help.\n'
           'Usage:\n'
           '\tpython tango_ass_2_taxa_freq.py -n nodes_file -d division -i tango_output\n'
           '\t'
    )


NODESFILE = ""
basename = ""
tango_output = ""
try:
    opts, args = getopt.getopt(sys.argv[1:], "hd:n:i:h")
except getopt.GetoptError, err:
    print str(err)
    usage()
    sys.exit()
for o, a in opts:
    if o == "-h":
        usage()
        sys.exit()
    elif o == "-n":
        NODESFILE = a
    elif o == "-d":
        basename = a
    elif o == "-i":
        tango_output = a+"_0.5"
    else:
        assert False, "Unhandled option."

if NODESFILE != "":
    print NODESFILE
else:
    print "No dump file"
    exit()

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
for line in open(NODESFILE):
    line = line.strip()
    fields = map(strip, line.split("|"))
    s = fields[0].split("#")
    nodeid = s[0]
    c = fields[1].split("#")
    parent = c[0]
    node2name[s[0]] = s[1]
    node2order[s[0]] = s[2]
    node2parent[nodeid] = parent


##################################
# ANALISI DEI DATI			     #
##################################
#analisi risultati tango
out_acc = set()
seq2taxa = {}
result = open(tango_output)
for line in result.readlines():
    line = line.strip()
    s = line.split("\t")
    out_acc.add(s[0])
    path = s[2].split(";")
    if node2order[path[0]] == "GB acc":
        seq2taxa[s[0]] = node2parent[path[0]]
    else:
        seq2taxa[s[0]] = path[0]
result.close()


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


