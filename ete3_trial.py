from ete3 import Tree, TreeStyle, faces, AttrFace, NodeStyle
from ete3 import NCBITaxa
from ete3 import nexml
from ete3 import Phyloxml, phyloxml
from Bio import Phylo
from ete3 import PhyloTree
import os

'''
items = [['a', 'b', 'c', 'd'],['a', 'b', 'c'], ['a','b','c','d','e'], ['a','b','c','d','e', 'f'], ['a','b','c','d','e', 'g'], ['a','b','c','d','h']]
'''
with open('output.txt',encoding = 'utf-8') as fp:
                                    os.mkdir("trees")
                                    text = fp.readlines()
                                    topologyfeeder = []
                                    for i in text:
                                          id = i[6:i.find(']')]
                                          id = id.replace("'", "")
                                          topologyfeeder.append(str(id))
                                    print(topologyfeeder)
                                  
                                    ncbi = NCBITaxa()
                                    
                                    tree = ncbi.get_topology(topologyfeeder, intermediate_nodes=True)
                                    #print(tree.get_ascii(attributes=["sci_name", "rank"]))
                                    
                                    
                                   
                                    ts = TreeStyle()
                                    ts.show_leaf_name = False
                                    ts.mode = "c"
                                    ts.root_opening_factor = 1
                                    ts.arc_start = -180 # 0 degrees = 3 o'clock
                                    ts.arc_span = 360
                                    tree.show(tree_style=ts)
                                    tree.write(format = 1, outfile = "trees/new_tree.nwk")
                                    tree2 = PhyloTree("trees/new_tree.nwk", sp_naming_function=lambda name: name)
                                    tax2names, tax2lineages, tax2rank = tree2.annotate_ncbi_taxa()
                                    tree2.write(format = 1, outfile = "trees/new_tree2.nwk")
                                    
tree = Phylo.read("trees/new_tree.nwk", "newick")

Phylo.convert("trees/new_tree.nwk", "newick", "trees/new_tree.xml", "nexml")
