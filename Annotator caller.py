#Annotator caller
import PySimpleGUI as sg
from edited_pipeline import Annotator as ann
'''
jsondic = sg.popup_get_file("Location of NCBI_TAX_DIC.json:")
sg.popup('You entered', jsondic)
inputfiles = sg.popup_get_folder("Folder of input bioc files:")
sg.popup("You entered", inputfiles)
outputfiles = sg.popup_get_folder("Directory to contain output json files:")
sg.popup("You entered", outputfiles)
result = ann(jsondic, inputfiles, outputfiles)
result.initialsteps()
'''
result = ann('New_Pipeline_2023/newdic.json', 'output', 'New_Pipeline_2023')
result.initialsteps()