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
                                        skipper = False
                                        for index, word in enumerate(wordlist):
                                            wordaslist = list(word)
                                            finalword = []
                                            for i in wordaslist:
                                                 if i == '.' or i == ',':
                                                      continue
                                                 else:
                                                      finalword.append(i)
                                            finalword = "".join(finalword)
                                            wordlist[index] = finalword
                                        print(wordlist)
                                        for index, word in enumerate(wordlist):
                                            if skipper == True:
                                                skipper = False
                                                continue
                                            
                                            else:
                                                if word in CleanNames:
                                                        #Exact_match
                                                        match = next((l for l in dict_data if l['CleanName'] == word), None)
                                                        if index < len(wordlist) -1:
                                                            if match['TaxRank'] == "genus" and (wordlist[index + 1] not in ['sp', 'spp', 'genus', 'gen']):
                                                                possible_species = word + " " + str(wordlist[index + 1])
                                                                if possible_species in CleanNames:
                                                                    match = next((l for l in dict_data if l['CleanName'] == possible_species), None)
                                                                    modifier = "species"
                                                                    self.AddAnnotation(match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext)
                                                                    skipper = True
                                                                else:
                                                                    self.AddAnnotation(match, self.count, m, " ", taxa_per_file, sentenceoffset, offsetoftext)
                                                                    skipper = False
                                                                    
                                                            elif wordlist[index + 1] in ['sp', 'spp']:
                                                                modifier = "species"
                                                                self.AddAnnotation(match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext)
                                                                skipper = False
                                                            elif wordlist[index + 1] in ['genus', 'gen']:
                                                                modifier = "genus"
                                                                self.AddAnnotation(match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext)
                                                                skipper = False
                                                            else:
                                                                 self.AddAnnotation(match, self.count, m, " ", taxa_per_file, sentenceoffset, offsetoftext)
                                                                 skipper = False
                                                        else:
                                                             self.AddAnnotation(match, self.count, m, " ", taxa_per_file, sentenceoffset, offsetoftext)
                                                             skipper = False
     
                                                     
                                                else:   
                                                        possible = []
                                                        for cn in CleanNames:
                                                            cn.split(" ")
                                                            if word in cn:
                                                                #Cleannames may be multiple words. We do not want words within words. eg: tar in starship. 
                                                                    possible.append(cn)
                                                        #Possible is a list of clean names which contain a word from the text.  If there is a possible list, do the following:     
                                                              
                                                        if len(possible) == 0:
                                                             #Try latin plural endings
                                                             if word.endswith('ae'):
                                                                newword = re.sub('ae$', 'a', word)
                                                             elif word.endswith('i'):
                                                                newword = re.sub('i$', 'us', word)
                                                             elif word.endswith('a'):
                                                                newword= re.sub('a$', 'um', word)
                                                             else:
                                                                  continue
                                                             if newword in CleanNames:
                                                                  match = next((l for l in dict_data if l['CleanName'] == newword), None)
                                                                  self.AddAnnotation(match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext)
                                                                  skipper = False
                                                                          
                                                             for cn in CleanNames:
                                                                    cn.split(" ")
                                                                    if newword in cn:
                                                                         possible.append(cn)
                                                        else: continue     
                                                                     
                                                        if len(possible) >=1:
                                                            for p in possible:
                                                                p= p.split(" ")
                                                                longest = max(0, len(p))
                                                                 # Scans for words nearby a recognised word. Accounts for multipart words
                                                            scanstart = index - longest
                                                            scanend = index + longest
                                                            for n in range(scanstart, scanend):
                                                                section = wordlist[n: n + longest:1]
                                                                section = " ".join(section)
                                                                
                                                                if section in possible:

                                                                        
                                                                        for d in dict_data:
                                                                                if d['CleanName'] == section:
                                                                                    match = d
                                                                                    self.AddAnnotation(match, self.count, m, " ", taxa_per_file, sentenceoffset, offsetoftext)
                                                                                    skipper = False
                                                                                    break
                                                                                else:
                                                                                    continue
                                                                else: # Format E. coli
                                                                     section = section.split(" ")
                                                                     if (len(section) ==2) :
                                                                          if section[0][0].isupper()  and (section[0][1] == '.'):
                                                                            for p in possible:
                                                                                if section[0][0] == p[0] and p in taxa_per_file:
                                                                                     match = next((l for l in dict_data if l['CleanName'] == p), None)
                                                                                     modifier = "species"
                                                                                     self.AddAnnotation(match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext)
                                                                                     skipper = False
                                                                          
                                                                          
                                                            
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

        def AddAnnotation(self, match, count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext):
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
                    taxa_per_file.append(str(match['CleanName']))
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
                    taxa_per_file.append(match['CleanName'])
    
            m['annotations'].append(dictannot)


             

# Latin noun endings
'''
G1      F       puella      puellae
G2      M       servus      servi
G2      N       templum     templa
G3      M/F     rex         reges
G3      N       corpus      corpora
G3 are difficult because the noun stem changes. Most of the focus should be therefore on G1 and G2. 
G4 and G5 are uncommon and don't change much in nominative singular and plural. 
'''
