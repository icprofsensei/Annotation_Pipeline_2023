from ete3 import Tree, TreeStyle
from ete3 import NCBITaxa
import sys
from ete3 import nexml


'''
items = [['a', 'b', 'c', 'd'],['a', 'b', 'c'], ['a','b','c','d','e'], ['a','b','c','d','e', 'f'], ['a','b','c','d','e', 'g'], ['a','b','c','d','h']]
'''
with open('output.txt',encoding = 'utf-8') as fp:
                                    text = fp.readlines()
                                    topologyfeeder = []
                                    for i in text:
                                          id = i[6:i.find(']')]
                                          id = id.replace("'", "")
                                          topologyfeeder.append(str(id))
                                    print(topologyfeeder)
                                    ncbi = NCBITaxa()
                                    
                                    tree = ncbi.get_topology(topologyfeeder, intermediate_nodes=True)
                                    print(tree.get_ascii(attributes=["sci_name", "rank"]))
                                    tree.write(format = 100, outfile = "new_tree.nw")
                                    
                                    ts = TreeStyle()
                                    ts.show_leaf_name = True
                                    ts.mode = "c"
                                    ts.arc_start = -180 # 0 degrees = 3 o'clock
                                    ts.arc_span = 180
                                    tree.show(tree_style=ts)

# Start with the shortest sequence = Basis
# Find trees with length = shortest + 1
# Add their last items inside brackets 