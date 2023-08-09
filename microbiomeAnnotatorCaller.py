#Annotator caller
import PySimpleGUI as sg
from microbiomeAnnotator import Annotator as ann
'''
jsondic = sg.popup_get_file("Location of NCBI_Tax_dictionary.json:")
sg.popup('You entered', jsondic)
inputfiles = sg.popup_get_folder("Folder of input bioc files:")
sg.popup("You entered", inputfiles)
outputfiles = sg.popup_get_folder("Directory to contain output json files:")
sg.popup("You entered", outputfiles)
result = ann(jsondic, inputfiles, outputfiles, 0)
result.initialsteps()
'''
result = ann('NCBI_tax_dictionary8.json', 'Testset3', 'results', 0)
result.initialsteps()
