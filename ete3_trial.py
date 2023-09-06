from ete3 import Tree, TreeStyle, faces, AttrFace, NodeStyle
from ete3 import NCBITaxa
from ete3 import nexml
from ete3 import Phyloxml, phyloxml
from Bio import Phylo
from ete3 import PhyloTree
import os
import requests
import math

'''
items = [['a', 'b', 'c', 'd'],['a', 'b', 'c'], ['a','b','c','d','e'], ['a','b','c','d','e', 'f'], ['a','b','c','d','e', 'g'], ['a','b','c','d','h']]
'''

class Tree:
        def __init__(self, items_to_find, directorypath):
            #Initialise inputs
            self.items_to_find = items_to_find
            self.directorypath = directorypath
        def listmaker(self, listtobeprocessed, allitems):
                ncbi = NCBITaxa()
                for itf in listtobeprocessed:
                        for i in ncbi.get_lineage(itf):
                                allitems.append(str(i))
                                #print(ncbi.get_lineage(itf))
                allitems = list(dict.fromkeys(allitems))
                return allitems
        def colourselecter(self):
                allitems = self.listmaker(self.items_to_find, [])
                iddict=dict.fromkeys(allitems,0)
                #print(iddict)
                ncbi = NCBITaxa()
                with open ('childnodes.txt',encoding = 'utf-8') as cn:
                        text = cn.readlines()
                        childnodedict = dict()
                        for i in text:
                                item = i.split(" ")
                                if item[1] == "\n":
                                        item[1] = 1
                                else:
                                        item[1] = item[1].strip("\n")
                                key = item[0]
                                value = item[1]
                                childnodedict[key] = value
                #print(childnodedict)
                for itf in self.items_to_find:
                        reversedls = ncbi.get_lineage(itf)[::-1]
                        #print(reversedls)
                        factor = 1
                        for index, i in enumerate(reversedls):
                                        #print(i)
                                        if index == 0:
                                                #print("Factor:", 1)
                                                iddict[str(i)] += 1 
                                        else:
                                                if str(i) in childnodedict.keys():

                                                        newfactor = 1/ int(childnodedict[str(i)])
                                                        factor = factor * newfactor
                                                        #print("Factor:", factor)
                                                        iddict[str(i)] += factor
                                                else:
                                                        iddict[str(i)] += factor 
                total = max(iddict.values())
                viridis = ['#fde725',
'#f8e621',
'#f1e51d',
'#ece51b',
'#e5e419',
'#dfe318',
'#d8e219',
'#d0e11c',
'#cae11f',
'#c2df23',
'#bddf26',
'#b5de2b',
'#addc30',
'#a8db34',
'#a0da39',
'#9bd93c',
'#93d741',
'#8ed645',
'#86d549',
'#7fd34e',
'#7ad151',
'#73d056',
'#6ece58',
'#67cc5c',
'#60ca60',
'#5cc863',
'#56c667',
'#52c569',
'#4cc26c',
'#48c16e',
'#42be71',
'#3dbc74',
'#3aba76',
'#35b779',
'#32b67a',
'#2eb37c',
'#2ab07f',
'#28ae80',
'#25ac82',
'#24aa83',
'#22a785',
'#20a486',
'#1fa287',
'#1fa088',
'#1f9e89',
'#1e9b8a',
'#1f998a',
'#1f968b',
'#20938c',
'#20928c',
'#218f8d',
'#228d8d',
'#238a8d',
'#24878e',
'#25858e',
'#26828e',
'#26818e',
'#277e8e',
'#287c8e',
'#29798e',
'#2a768e',
'#2b748e',
'#2c718e',
'#2d708e',
'#2e6d8e',
'#306a8e',
'#31688e',
'#32658e',
'#33638d',
'#34608d',
'#365d8d',
'#375b8d',
'#38588c',
'#39558c',
'#3b528b',
'#3c508b',
'#3d4d8a',
'#3e4989',
'#3f4788',
'#414487',
'#424186',
'#433e85',
'#443a83',
'#453882',
'#463480',
'#46327e',
'#472e7c',
'#472c7a',
'#482878',
'#482475',
'#482173',
'#481d6f',
'#481b6d',
'#481769',
'#471365',
'#471063',
'#460b5e',
'#46085c',
'#450457',
'#440154']
                reverseviridis = viridis[::-1]
                for key, value in iddict.items():
                        index = (value / total) * 100
                        index = math.ceil(index)
                        iddict[key] = reverseviridis[index - 1]
                        
                return(iddict) 
                                        

        def layoutfunc(self, node):
                  node.complete_branch_lines_when_necessary = False
                  node.optimal_scale_level = "full"
                  node.guiding_lines_type = 0
                  node.extra_branch_line_type = 0
                  mydict = self.colourselecter()
                  
                  node.img_style["hz_line_type"] = 0
                  if node.get_children() == [] or node.name not in mydict.keys():
                          node.img_style["hz_line_color"] = "#ffffff"
                  nohorline = False
                  
                  
                  if node.name in mydict.keys():
                        #node.img_style["size"] = 0.1
                        #node.img_style["shape"] = "circle"
                        node.img_style["fgcolor"] = mydict[node.name]
                        node.img_style["vt_line_color"] = mydict[node.name]
                        node.img_style["hz_line_color"] = mydict[node.name]
                        #faces.add_face_to_node(AttrFace("name", fsize = 25), node, column=0)
                        if node.get_children == []:
                                for i in node.get_children():
                                        if i in mydict.keys():
                                                continue
                                        else:
                                                nohorline == True
        
                  else: 
                        
                        node.img_style["size"] = 0
                        node.img_style["fgcolor"] = "#ffffff"
                        node.img_style["bgcolor"] = "#ffffff"
                        node.support = 0
                        node.distance = 0
                        node.img_style["vt_line_color"] = "#ffffff"
                        node.img_style["hz_line_color"] = "#ffffff"
                        
                        if node.get_children() != []:
                                nohorline == True
                                        
                  if nohorline == True:
                              node.img_style["hz_line_color"] = "#ffffff"      
                               
        def Maker(self):
                
                  with open('output.txt',encoding = 'utf-8') as fp:
                                                      os.mkdir(self.directorypath + "trees")
                                                      text = fp.readlines()
                                                      topologyfeeder = []
                                                      for i in text:
                                                            id = i[6:i.find(']')]
                                                            id = id.replace("'", "")
                                                            topologyfeeder.append(str(id))
                                                      #print(topologyfeeder)
                                                
                                                      ncbi = NCBITaxa()
                                                      
                                                      tree = ncbi.get_topology(topologyfeeder, intermediate_nodes=True)
                                                      
                                                      #print(tree.get_ascii(attributes=["sci_name", "rank"]))
                                                      ts = TreeStyle()
                                                      ts.show_leaf_name = False
                                                      ts.mode = "c"
                                                      ts.root_opening_factor = 0
                                                      ts.arc_start = -180 # 0 degrees = 3 o'clock
                                                      ts.arc_span = 360
                                                      #ts.branch_vertical_margin = 0.5
                                                      
                                                      #ts.min_leaf_separation = 6
                                                      ts.layout_fn = self.layoutfunc
                                                      tree.show(tree_style=ts)
                                                      tree.write(format = 1, outfile = self.directorypath + "trees/new_tree.nwk")
                                                      tree.render(self.directorypath + "trees/img_tree.svg", w= 3600, units = 'px', tree_style = ts)
                                                      Phylo.convert(self.directorypath + "trees/new_tree.nwk", "newick", self.directorypath + "trees/new_tree.xml", "nexml")
                                                      '''
                                                      tree2 = PhyloTree("trees/new_tree.nwk", sp_naming_function=lambda name: name)
                                                      tax2names, tax2lineages, tax2rank = tree2.annotate_ncbi_taxa()
                                                      tree2.write(format = 1, outfile = "trees/new_tree2.nwk")
                                                      '''


