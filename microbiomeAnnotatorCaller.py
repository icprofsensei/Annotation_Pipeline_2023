#Annotator caller
import PySimpleGUI as sg
from microbiomeAnnotator_condensed import Annotator as ann


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
tree = sg.popup_get_file("Location of childnodes.txt file (if answered YES previously, else click OK):")
sg.popup('You entered', tree)
result = ann(jsondic, inputfiles, outputfiles, 0, values, values2, tree)
result.initialsteps()

'''

result = ann('NCBI_tax_dictionary8.json', 'testset', 'results', 0, 'ALL', 'YES', 'childnodes.txt')
result.initialsteps()
