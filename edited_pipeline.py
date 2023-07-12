import numpy as np
import pandas as pd
import json
import re
import time
from collections import Counter
import os
import datetime
import os.path
from alive_progress import alive_bar


class Annotator:
        def __init__(self, dic_directory, input_directory, output_directory):
            #Initialise inputs
            self.dic_directory = dic_directory          #dic.json
            self.input_directory = input_directory      #Input bioc files
            self.output_directory = output_directory    #Output file directory where Annotated_ouput is created. 

        def initialsteps(self):
            os.mkdir(self.output_directory + "/Annotated_output2")
            f1 = open(self.dic_directory, encoding = 'utf-8') 
            if os.path.isfile(self.dic_directory) == True:
                 print("Dictionary loaded")
            dict_data = json.load(f1)
            CleanNames = []
            for i in dict_data:
                CleanNames.append(i['CleanName'])
            CleanNames =list(dict.fromkeys(CleanNames))
            print(CleanNames)
            list1=[]
            start_time = datetime.datetime.now() #start time

            message = '' #message with start and stop time

            all_files = os.listdir(self.input_directory) 

            PMC_files=[]  
            documents,passages=[],[]
            id=0

#PMC_files contains all files in the input which begin with PMC. 
            for n in all_files:   
                if n.startswith('PMC') and n.endswith('bioc.json'):  
                    PMC_files.append(n)


#For each PMC file, the data is loaded as a json and my_list is made. Text is under documents --> passages --> annotations --> text --> word     
            
                                      
            for in_file in PMC_files: 
                    with open(self.input_directory+ "/" + in_file , encoding = 'utf-8') as m_file:
                        if os.path.isfile(self.input_directory+ "/" + in_file) == True:
                                print("Input file found")
                        data = json.load(m_file)
                        my_list=[]
                        documents=data['documents'] 
                        for j in documents:
                                passages= j['passages'] 
                                total = len(passages)
                                with alive_bar(total) as bar: 
                                    for m in passages:
                                        m['annotations']=[]      
                                        textsection=m['text']
                                        wordlist=textsection.split(" ")
                                        #Iterates over a list of words in the wordlist, taken from the text section.
                                        
                                        for word in wordlist:
                                            #if not word.startswith(" "):
                                                #if word[0].isupper() == True:
                                                        if word in CleanNames:
                                                            match = next((l for l in dict_data if l['CleanName'] == word), None)
                                                            my_list.append(word, {match['CleanName'], match['TaxRank']})
                                                            m['annotations'].append({
                                                                                        "text":textsection,
                                                                                        "infons":{
                                                                                            "name":word,
                                                                                            "identifier": match['CleanName'] ,
                                                                                            "taxonomic rank": match['TaxRank'] ,
                                                                                            "annotator":"dhylan.patel21@imperial.ac.uk",
                                                                                            "date": time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()) ,
                                                                                        },
                                                                                        "id": match['TaxID'],
                                                                                        "locations":{
                                                                                            "length": len(textsection),
                                                                                            "parent taxonomic id": match['ParentTaxID'],
                                                                                            
                                                                                        }
                                                                                    })
                                                        else:
                                                            continue
                                                    
                                        bar()
                                                
                                        
                #saving it as a json file 
                                print(my_list)
                                list1.append(my_list)
                                for i in data: 
                                    a_file = open(self.output_directory + "/Annotated_output2"+ "/" +str(in_file), "w")
                                    json.dump(data, a_file, indent = 4)
                                    a_file.close()
             
            stop_time = datetime.datetime.now() #stop time
            print(id)
            message = 'Start time is ' + str(start_time) + '\n' + 'Stop time is ' + str(stop_time)  

            with open(self.output_directory + "/Annotated_output2"+ "/" + 'runtime11.log', 'a') as time_file:   #start and stop time are written into a runtime.log file
                
            
                time_file.write(message)

            with open(self.output_directory + "/Annotated_output2"+ "/" + 'runtime12new.log', 'a') as time_file:   #start and stop time are written into a runtime.log file
                
                time_file.write("This file contains all the kingdoms and taxonomic ranks for each species found.")
                time_file.write(str(list1))

