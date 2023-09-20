# Caller
from ete3_trial import TreeMaker as T
call = T(['43675', '644652'], 'results2/', 'childnodes.txt', 'my choice')
call.Maker()
'''
from folder_parser import species_miner as sm
result = sm('All_annotated_files', 'introduction')
result.Miner()
result.Analyser()
result.Analyser()

from childnodes import Childnodes as C
call = C('childnodes.txt', ['9606', '9771'], "f55726c2c32772c2b82304814b30148aff07")
call.filerefresh()'''