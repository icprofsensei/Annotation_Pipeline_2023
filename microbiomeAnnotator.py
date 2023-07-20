import numpy as np
import pandas as pd
import json
import re
import time
import os
import datetime
import os.path
from alive_progress import alive_bar


class Annotator:
        def __init__(self, dic_directory, input_directory, output_directory, count):
            #Initialise inputs
            self.dic_directory = dic_directory          #newdic.json
            self.input_directory = input_directory      #Input bioc files
            self.output_directory = output_directory    #Output file directory where Annotated_ouput is created. 
            self.count = count

        def initialsteps(self):
            folder = '/Annotated_output_' + str(time.strftime('%Y-%m-%d_%H-%M-%S',time.localtime()))
            os.mkdir(self.output_directory + folder)
            opendic = open(self.dic_directory, encoding = 'utf-8') 
            if os.path.isfile(self.dic_directory) == True:
                 print("Dictionary loaded")
            dict_data = json.load(opendic)
            CleanNames = []
            for i in dict_data:
                CleanNames.append(i['CleanName'])
            CleanNames =list(dict.fromkeys(CleanNames))
            taxalist=[]
            start_time = datetime.datetime.now() #start time
            message = '' #message with start and stop time
            all_files = os.listdir(self.input_directory) 

#PMC_files contains all files in the input which begin with PMC. 
            PMC_files=[]  
            for n in all_files:   
                if n.startswith('PMC') and n.endswith('bioc.json'):  
                    PMC_files.append(n)


#For each PMC file, the data is loaded as a json and my_list is made. Text is under documents --> passages --> annotations --> text --> word     
            
            for i in range(len(PMC_files)):
                    in_file = PMC_files[i]
                    print("Annotating file: ", i + 1, "of: ", len(PMC_files))                         
                    count = str(self.count)
                    with open(self.input_directory+ "/" + in_file , encoding = 'utf-8') as m_file:
                        if os.path.isfile(self.input_directory+ "/" + in_file) == True:
                                print("Input file found")
                        data = json.load(m_file)
                        taxa_per_file=[]
                        documents = []
                        documents=data['documents'] 
                        for j in documents:
                                passages = []
                                passages= j['passages'] 
                                total = len(passages)
                                with alive_bar(total) as bar: 
                                    for m in passages:
                                        m['annotations']=[]      
                                        textsection=m['text']
                                        offsetoftext = m['offset']
                                        wordlist=textsection.split(" ")
                                    #Iterates over a list of words in the wordlist, taken from the text section.
                                        sentenceoffset = 0
                                        for index, word in enumerate(wordlist):
                                            for i in word:
                                                if i.isalnum() == True:
                                                    if i.isupper():
                                                        break
                                                    else:
                                                        word = word.rstrip(".,")
                                                        wordlist[index] = word.rstrip(".,")
                                                else: 
                                                    continue
                                        for index, word in enumerate(wordlist):
                                            if len(word) >= 4:
                                                if word in CleanNames: 
                                                        #Exact_match
                                                        match = next((l for l in dict_data if l['CleanName'] == word), None)
                                                        if match['TaxRank'] == "genus":
                                                             possible_species = word + " " + str(wordlist[index + 1])
                                                             print(possible_species)
                                                             if possible_species in CleanNames:
                                                                  print("found")
                                                                  continue
                                                             else:
                                                                  self.AddAnnotation(match, self.count, index, m, " ", taxa_per_file, sentenceoffset, offsetoftext)
                                                        elif wordlist[index + 1].rstrip(".,") in ['[sp]', '[spp]', '[species]', 'sp', 'spp']:
                                                            modifier = "species"
                                                            self.AddAnnotation(match, self.count, index, m, modifier, taxa_per_file, sentenceoffset, offsetoftext)
                                                        elif wordlist[index + 1].rstrip(".,") in ['[gen]', '[genus]', 'genus', 'gen']:
                                                            modifier = "genus"
                                                            self.AddAnnotation(match, self.count, index, m, modifier, taxa_per_file, sentenceoffset, offsetoftext)
                                                        else:
                                                             self.AddAnnotation(match, self.count, index, m, " ", taxa_per_file, sentenceoffset, offsetoftext)
                                        
                                                else:   
                                                        possible = []
                                                        for cn in CleanNames:
                                                            if word in cn:
                                                                x = cn.find(word)
                                                                y = len(word)
                                                                if  (cn.startswith(word) and cn[x+y] == " ") or (cn.endswith(word) and cn[x-1] == " ") or (cn[x-1] == " " and cn[x+y] ==" ") : #Cleannames may be multiple words. We do not want words within words. eg: tar in starship. 
                                                                    possible.append(cn)
                                                        #Possible is a list of clean names which contain a word from the text.  If there is a possible list, do the following:             
                                                        if len(possible) >=1:
                                                            for p in possible:
                                                                p= p.split(" ")
                                                                longest = max(0, len(p))
                                                                
                                                            scanstart = index - longest
                                                            scanend = index + longest
                                                            for n in range(scanstart, scanend):
                                                                section = wordlist[n: n + longest:1]
                                                                section = " ".join(section)
                                                                #print(section)
                                                                if section in possible:

                                                                        #print("Identified", section)
                                                                        for d in dict_data:
                                                                                if d['CleanName'] == section:
                                                                                    match = d
                                                                                    self.AddAnnotation(match, self.count, index, m, " ", taxa_per_file, sentenceoffset, offsetoftext)
                                                                        break
                                                                else:
                                                                            continue
                                            sentenceoffset += (len(word)+ 1)
                                        bar()

                                                          

                                taxa_per_file = {*taxa_per_file}
                                taxa_per_file = list(taxa_per_file)
                                for j in taxa_per_file:
                                      print(j)
                                taxalist.append(taxa_per_file)
                                for i in data: 
                                    a_file = open(self.output_directory + folder + "/" +str(in_file), "w")
                                    json.dump(data, a_file, indent = 4)
                                    a_file.close()
             
            stop_time = datetime.datetime.now() #stop time
            message = 'Start time is ' + str(start_time) + '\n' + 'Stop time is ' + str(stop_time)  

            with open(self.output_directory + folder + "/" + 'runtime.log', 'a') as time_file:   #start and stop time are written into a runtime.log file
                
            
                time_file.write(message)

            with open(self.output_directory + folder + "/" + 'taxafound.log', 'a') as time_file:   #All taxa found are added to taxafound.log file.
                
                time_file.write("This file contains all the kingdoms and taxonomic ranks for each species found.")
                time_file.write(str(taxalist))

        def AddAnnotation(self, match, count, index, m, modifier, taxa_per_file, sentenceoffset, offsetoftext):
            self.count = int(count) + 1
            
            if modifier != " ": 
                    dictannot = {
                                        "text":match['CleanName'],
                                        "infons":{
                                            "identifier":  match['TaxID'] ,
                                            "type": modifier ,
                                            "annotator":"dhylan.patel21@imperial.ac.uk",
                                            "date": time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()) ,
                                            "parent_taxonomic_id": "Modified type - not identified"
                                        },
                                        "id": count,
                                        "locations":{
                                            "length": len(match['CleanName']),
                                            "offset": sentenceoffset + offsetoftext ,
                                            
                                        }
                                    }
                    taxa_per_file.append(str(match['CleanName']) + " " + str(match['CleanName']) + " " + modifier)
            else:
                    dictannot = {
                                        "text":match['CleanName'],
                                        "infons":{
                                            "identifier": match['TaxID'] ,
                                            "type": match['TaxRank'] ,
                                            "annotator":"dhylan.patel21@imperial.ac.uk",
                                            "date": time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()) ,
                                            "parent_taxonomic_id": match['ParentTaxID']
                                        },
                                        "id": count,
                                        "locations":{
                                            "length": len(match['CleanName']),
                                            "offset": sentenceoffset + offsetoftext ,
                                            
                                        }
                                    }
                    taxa_per_file.append(match['CleanName'] + " " + str(match['CleanName']) + " " + str(match['TaxRank']))
    
            m['annotations'].append(dictannot)