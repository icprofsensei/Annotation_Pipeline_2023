import numpy as np
import pandas as pd
import json
import re
import time
import os
import datetime
import os.path

from collections import Counter
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
            cncounts = dict(Counter(CleanNames))
            strains = {key:value for key, value in cncounts.items() if value > 1}
            CleanNames =list(dict.fromkeys(CleanNames))
            taxalist=[]
            start_time = datetime.datetime.now() 
            message = '' 
            all_files = os.listdir(self.input_directory) 
            PMC_files=[]  
            for n in all_files:   
                if  n.endswith('bioc.json'):  
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
                        
                        documents = []
                        documents=data['documents'] 
                        for j in documents:
                                passages = []
                                passages= j['passages'] 
                                total = len(passages)
                                with alive_bar(total) as bar: 
                                    taxa_per_file=[]
                                    for m in passages:
                                        m['annotations']=[]      
                                        textsection=m['text']
                                        offsetoftext = m['offset']
                                        oglist = textsection.split(" ")
                                        wordlist=textsection.split(" ")
                                        needs_processing  = []
                                        sentenceoffset = 0
                                        skipper = False
                                        for index, word in enumerate(wordlist):
                                            wordaslist = list(word)
                                            finalword = []
                                            if len(wordaslist)!=2:   
                                                for i in wordaslist:
                                                    if i == '.' or i == ',':
                                                        continue
                                                    else:
                                                        finalword.append(i)
                                                finalword = "".join(finalword)
                                                wordlist[index] = finalword
                                            else:
                                                 if wordaslist[0].isupper() and wordaslist[1] == '.':
                                                      continue
                                                 
                                        for index, word in enumerate(wordlist):
                                            if len(wordlist[index]) < len(oglist[index]): # This adjusts the offset based on removing '.' and ',' from words. Won't work if someone writes  .  or , on its own
                                                 sentenceoffset += 1
                                            idinuse = []
                                            newword = ""
                                            newword = self.CheckLatin(word, newword)
                                            duptxids = []
                                            

                                            if skipper == True:
                                                sentenceoffset += (len(word)+ 1) 
                                                skipper = False
                                                continue
                                            
                                            else:
                                                if word in CleanNames:
                                                        # EXACT MATCH
                                                        match = next((l for l in dict_data if l['CleanName'] == word), None)
                                                        
                                                        if index < len(wordlist) -1:
                                                            nextword = wordlist[index + 1]
                                                            newnextword = ""
                                                            newnextword = self.CheckLatin(nextword, newnextword) 
                                                            if match['TaxRank'] == "genus" and (nextword not in ['sp', 'spp', 'genus', 'gen']):
                                                                possible_species = word + " " + str(nextword)
                                                                possible_plural = word + " " + str(newnextword)
                                                                #Genus species - Do not annotate the next word 
                                                                if possible_species in CleanNames:
                                                                    match = next((l for l in dict_data if l['CleanName'] == possible_species), None)
                                                                    modifier = "species"
                                                                    self.AddAnnotation(match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing)
                                                                    skipper = True
                                                                #Genus species(pl) - Do not annotate the next word
                                                                elif possible_plural in CleanNames:
                                                                     match = next((l for l in dict_data if l['CleanName'] == possible_plural), None)
                                                                     modifier = "species"
                                                                     self.AddAnnotation(match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing)
                                                                     skipper = True
                                                                # Genus - Only one word, so continue to annotate the next word
                                                                else:
                                                                    self.AddAnnotation(match, self.count, m, " ", taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing)
                                                                    skipper = False   
                                                            # [Any] sp - Do not annotate the next word
                                                            elif nextword in ['sp', 'spp']:
                                                                modifier = "species"
                                                                self.AddAnnotation(match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing)
                                                                skipper = True
                                                            # [Any] genus - Do not annotate the next word
                                                            elif nextword in ['genus', 'gen']:
                                                                modifier = "genus"
                                                                self.AddAnnotation(match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing)
                                                                skipper = True
                                                        # [Any] Only one word, so continue to the next word.   
                                                        else:
                                                             self.AddAnnotation(match, self.count, m, " ", taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing)
                                                             skipper = False
                                                #  [Any pl] Only one word, so continue to annotate the next word - EXACT MATCH
                                                elif newword in CleanNames: 
                                                     match = next((l for l in dict_data if l['CleanName'] == newword), None)
                                                     self.AddAnnotation(match, self.count, m, " ", taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing)
                                                     skipper = False
                                                     
                                                else:   
                                                    if len(word) >= 4:
                                                        possible = []
                                                        #Checks if word matches a word within a clean name.
                                                        for cn in CleanNames:
                                                            cn = cn.split(" ")
                                                            if word in cn:
                                                                    possible.append(" ".join(cn))
                                                            elif newword in cn:
                                                                  possible.append(" ".join(cn))
                                                            else: continue     
                                                                     
                                                        if len(possible) >=1:
                                                            scanstart = index - 2
                                                            scanend = index + 2
                                                            for n in range(scanstart, scanend):
                                                                section = list(wordlist[n: n + 2:1])
                                                                if len(section) == 2:
                                                                    ls0 = len(section[0])
                                                                    # G. species
                                                                    if (ls0>0 and section[0][0].isupper()) and ((ls0==1) or (section[0][1] == '.')):
                                                                                for p in possible:
                                                                                    if section[0][0] == p[0] and p in taxa_per_file:
                                                                                        recurring = taxa_per_file.index(p)
                                                                                        match = next((l for l in dict_data if l['CleanName'] == taxa_per_file[recurring]), None)
                                                                                        modifier = "species"
                                                                                        self.AddAnnotation(match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing)
                                                                                        skipper = False
                                                                                    elif section[0][0] == p[0] and p not in taxa_per_file:
                                                                                        match = next((l for l in dict_data if l['CleanName'] == p), None)
                                                                                        modifier = "species"
                                                                                        if match['TaxID'] != idinuse:
                                                                                            self.AddAnnotation(match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing)
                                                                                            skipper = False
                                                                                        else:
                                                                                            continue
                                            sentenceoffset += (len(word)+ 1)  
                                               
                                        bar()
                                       # Post processing to remove ambiguity
                                        if needs_processing != []:
                                             #print("Needs post processing!")
                                             for ann, ian in enumerate(m['annotations']):
                                                if  '(' in ian['infons']['identifier']:
                                                     offset = ian['locations']['offset']
                                                     countog = ian["id"]
                                                     #
                                                     others = []
                                                     if ann == 0:
                                                        others.append(m['annotations'][ann + 1]['infons']['type'])
                                                     else:
                                                        if '(' not in m['annotations'][ann -1]['infons']['identifier']:
                                                            others.append(m['annotations'][ann - 1]['infons']['type'])
                                                        elif '(' not in m['annotations'][ann +1]['infons']['identifier']:
                                                             others.append(m['annotations'][ann + 1]['infons']['type'])
                                                        else: ian['infons']['type'] = 'unresolved'
                                                     loc_id = ian['infons']['identifier'].find('[')
                                                     unproc_str = ian['infons']['identifier'][loc_id:len(ian['infons']['identifier'])].split(" ")
                                                     possible_ids = []
                                                     items = []
                                                     for pid in unproc_str: 
                                                          pid = pid.strip("',[]")
                                                          possible_ids.append(pid)
                                                     for pid in possible_ids:
                                                          match = next((l for l in dict_data if l['TaxID'] == pid), None)
                                                          kingdic = {'NCBI:txid2': 'bacteria', 'NCBI:txid2157':'archaea', 'NCBI:txid4751':'fungi'}
                                                          kingdom = kingdic[match['KingdomID']]
                                                          identifier = kingdom + "_" + match['TaxRank'] 
                                                          item = [pid, identifier]
                                                          items.append(item)
                                        
                                                     for item in items:
                                                          if item[1] in others:
                                                               
                                                               match = next((l for l in dict_data if l['TaxID'] == item[0]), None)
                                                               
                                                               #print(match['CleanName'], item[1])
                                                               
                                                               dictannot = {
                                                                                "text":match['CleanName'],
                                                                                "infons":{
                                                                                    "identifier": item[0],
                                                                                    "type": item[1] ,
                                                                                    "annotator":"dhylan.patel21@imperial.ac.uk",
                                                                                    "date": time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()) ,
                                                                                    "parent_taxonomic_id": match['ParentTaxID']
                                                                                },
                                                                                "id": countog,
                                                                                "locations":{
                                                                                    "length": len(match['CleanName']),
                                                                                    "offset": offset ,
                                                                                    
                                                                                }
                                                                            }
                                                               m['annotations'].append(dictannot)            
                                                                
                                                          else:
                                                               continue
                                        m['annotations'].sort(key = lambda e: (e["id"]))    
                                        for ian in m['annotations']:
                                             if  '(' in ian['infons']['identifier'] and ian['infons']['type'] != 'unresolved':
                                                  m['annotations'].remove(ian)
                                             else:
                                                  continue
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

        def AddAnnotation(self, match, count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing):
            if match['TaxID'] not in idinuse:
                repeats= "", 
                self.count = int(count) + 1
                kingdic = {'NCBI:txid2': 'bacteria', 'NCBI:txid2157':'archaea', 'NCBI:txid4751':'fungi'}
                kingdom = kingdic[match['KingdomID']]
                if match['CleanName'] in strains:
                    repeats = int(strains[match['CleanName']])
        
                    for i in dict_data:
                        if i['CleanName'] == match['CleanName']:
                            duptxids.append(i['TaxID'])
                if repeats == "" or duptxids == []:
                     identifierstring = match['TaxID']
                else:
                     
                    identifierstring = match['TaxID'] + " (" + str(repeats) + ") " + str(duptxids)
                    needs_processing.append(identifierstring)
                
                if modifier != " ": 
                        dictannot = {
                                            "text":match['CleanName'],
                                            "infons":{
                                                "identifier": identifierstring,
                                                "type": kingdom + "_" + modifier ,
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
                        taxa_per_file.append(str(match['CleanName']))
                else:
                        dictannot = {
                                            "text":match['CleanName'],
                                            "infons":{
                                                "identifier": identifierstring,
                                                "type": kingdom + "_" + match['TaxRank'] ,
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
                idinuse.append(match['TaxID'])
                m['annotations'].append(dictannot)

        def CheckLatin (self, word, newword):
                if word.endswith('ae'):
                    newword = re.sub('ae$', 'a', word)
                elif word.endswith('i'):
                    newword = re.sub('i$', 'us', word)
                elif word.endswith('a'):
                    newword= re.sub('a$', 'um', word)
                return newword
  
             

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
