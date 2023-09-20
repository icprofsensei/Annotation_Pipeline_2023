#Annotator caller
import PySimpleGUI as sg
from microbiomeAnnotator_condensed import Annotator as ann

# (self, updatechildnodeyn, dic_directory, input_directory, output_directory, count, keyword, treeyn, treedir, ncbikey)
'''
jsondic = sg.popup_get_file("Location of NCBI_Tax_dictionary.json:")
sg.popup('You entered', jsondic)
inputfiles = sg.popup_get_folder("Folder of input bioc files:")
sg.popup("You entered", inputfiles)
outputfiles = sg.popup_get_folder("Directory to contain output json files:")
sg.popup("You entered", outputfiles)
layout = [[sg.Text('Please enter the section title being searched for. Type ALL for all section titles')],
[sg.Text('Section Title', size = (15,1)), sg.InputText()], [sg.Submit(), sg.Cancel()]]
window = sg.Window('Data Entry', layout)
event, values = window.read()
window.close()
layout2 = [[sg.Text('Type YES to produce a phylogenetic tree of entities.')],
[sg.Text('Section Title', size = (15,1)), sg.InputText()], [sg.Submit(), sg.Cancel()]]
window2 = sg.Window('Data Entry', layout2)
event, values2 = window2.read()
window2.close()
layout3 = [[sg.Text('Type YES to completely update the existing node descendent file - childnodes.txt (if answered YES previously, else click OK)')],
[sg.Text('Section Title', size = (15,1)), sg.InputText()], [sg.Submit(), sg.Cancel()]]
window3 = sg.Window('Data Entry', layout2)
event, values3 = window3.read()
window3.close()
cnodes = sg.popup_get_file("Location of existing childnodes.txt file")
sg.popup('You entered', cnodesdir)
layout4 = [[sg.Text('Enter the NCBI REST API key associated with your account')],
[sg.Text('Section Title', size = (15,1)), sg.InputText()], [sg.Submit(), sg.Cancel()]]
window4 = sg.Window('Data Entry', layout4)
event, values4 = window4.read()
window4.close()
result = ann(values3, jsondic, inputfiles, outputfiles, 0, values, values2, cnodesdir, values4)
result.initialsteps()

'''

result = ann('NO', 'NCBI_tax_dictionary8.json', 'testset', 'results', 0, 'ALL', 'YES', 'childnodes.txt', "f55726c2c32772c2b82304814b30148aff07")
result.initialsteps()
