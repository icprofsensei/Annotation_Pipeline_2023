#Annotator caller
import PySimpleGUI as sg
from microbiomeAnnotator_condensed import Annotator as ann
from tree_builder import TreeMaker as T
from childnodes import Childnodes as C

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
layoutcs = [[sg.Text('Type YES if you would like to turn case sensitivity ON')],
[sg.Text('Section Title', size = (15,1)), sg.InputText()], [sg.Submit(), sg.Cancel()]]
windowcs = sg.Window('Data Entry', layoutcs)
event, valuescs = windowcs.read()
window.close()
layout2 = [[sg.Text('Type YES to produce a phylogenetic tree of entities.')],
[sg.Text('Section Title', size = (15,1)), sg.InputText()], [sg.Submit(), sg.Cancel()]]
window2 = sg.Window('Data Entry', layout2)
event, values2 = window2.read()
window2.close()
layout3 = [[sg.Text('Type YES to completely update the existing node descendent file - childnodes.txt (if answered YES previously, else click OK)')],
[sg.Text('Section Title', size = (15,1)), sg.InputText()], [sg.Submit(), sg.Cancel()]]
window3 = sg.Window('Data Entry', layout3)
event, values3 = window3.read()
window3.close()
cnodesdir = sg.popup_get_file("Location of existing childnodes.txt file")
sg.popup('You entered', cnodesdir)
layout4 = [[sg.Text('Enter the NCBI REST API key associated with your account')],
[sg.Text('Section Title', size = (15,1)), sg.InputText()], [sg.Submit(), sg.Cancel()]]
window4 = sg.Window('Data Entry', layout4)
event, values4 = window4.read()
window4.close()
result = ann(jsondic, inputfiles, outputfiles, 0, values, valuescs)
result.initialsteps()
if type(values2) == dict:
                 
                values2 = str(values2[0])
if type(values3) == dict:
                 
                values3 =  str(values3[0])

if type(values4) == dict:
                 
                values4 = str(values4[0])

if values2 == 'YES':
        ids  = sg.popup_get_file("Location of ids.txt file")
        sg.popup('You entered', ids)
        if values3 == 'YES':
                print("updating the existing childnodes.txt file")
                result2 = C(cnodesdir, ids, values4)
                result2.filerefresh()
        print('Checking if any new species need to be added')
        result2 = C(cnodesdir, ids, values4)
        result2.updatenewspec()
        print('Check complete. Constructing trees.')
        result3 = T(ids, outputfiles, cnodesdir, inputfiles)
        result3.Maker()
'''
# ann(tax dictionary, input directory, output directory, count, section)
result = ann('NCBI_tax_dictionary8.json', 'testset', 'results', 0, 'ALL', 'NO')
result.initialsteps()
# C(childnodesdir, ids found text file, ncbiapikey)
#result2 = C('childnodes.txt', 'results/ALLAnnotated_output_2023-09-22_17-01-00/ids.txt', 'f55726c2c32772c2b82304814b30148aff07')
#result2.updatenewspec()
# T( ids found text file, output directory, childnodesdir, input directory)
# result3 = T('Testset4/ids.txt', 'results4', 'childnodes.txt', 'Testset4')
# result3.Maker()