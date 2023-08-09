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
                    print("Annotating file: ", i + 1, "of: ", len(PMC_files), in_file)    
                    base = i                   
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
                                    problemwords = []
                                    for m in passages:
                                        m['annotations']=[]      
                                        textsection=m['text']
                                        offsetoftext = m['offset']
                                        wordlist=textsection.split(" ")
                                        needs_processing  = []
                                        sentenceoffset = 0
                                        skipper = False
                                        
                                        for index, word in enumerate(wordlist):
                                            annot_stopper = False # Allows the use of multiple annotations for the same word. 
                                            idinuse = [] 
                                            duptxids = [] 
                                            if skipper == True:
                                                sentenceoffset += (len(word)+ 1) 
                                                skipper = False
                                                continue
                                            else:
                                                try:
                                                      
                                                    finalword = self.RemovePunc(word, [])             
                                                    wordlist[index] = finalword 
                                                    #print(finalword)
                                                    newword = ""
                                                    newword = self.CheckLatin(finalword, newword)
                                                    
                                                    possible = []
                                                    #Checks if word matches a word within a clean name.
                                                    for cn in CleanNames:
                                                        cn = cn.split(" ")
                                                        if finalword in cn:
                                                                possible.append(" ".join(cn))
                                                        else: continue  
                                                    if finalword in CleanNames:
                                                            
                                                            match = next((l for l in dict_data if l['CleanName'] == finalword), None)
                                                            
                                                            if index <= len(wordlist) -2:
                                                                            nextword = wordlist[index + 1]
                                                                            nextword = self.RemovePunc(nextword, [])  
                                                                            newnextword = ""
                                                                            newnextword = self.CheckLatin(nextword, newnextword) 
                                                                            possible_species = finalword + " " + str(nextword)
                                                                            possible_plural = finalword + " " + str(newnextword)
                                                                            
                                                                    
                                                                                        
                                                                            #Genus species - Do not annotate the next word 
                                                                            
                                                                            if possible_species in CleanNames:
                                                                                match = next((l for l in dict_data if l['CleanName'] == possible_species), None)
                                                                                modifier = "species"
                                                                                self.AddAnnotation(possible_species, match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing, base, annot_stopper)
                                                                                skipper = True
                                                                                annot_stopper = True
                                                                            #Genus species(pl) - Do not annotate the next word
                                                                            elif possible_plural in CleanNames:
                                                                                match = next((l for l in dict_data if l['CleanName'] == possible_plural), None)
                                                                                modifier = "species"
                                                                                self.AddAnnotation(possible_plural, match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing, base, annot_stopper)
                                                                                skipper = True
                                                                                annot_stopper = True
                                                                            # [Any] sp - Do not annotate the next word
                                                                            elif nextword in ['sp', 'spp']:
                                                                                modifier = "species"
                                                                                self.AddAnnotation(finalword, match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing, base, annot_stopper)
                                                                                skipper = True
                                                                                annot_stopper = True
                                                                            # [Any] genus - Do not annotate the next word
                                                                            elif nextword in ['genus', 'gen']:
                                                                                modifier = "genus"
                                                                                self.AddAnnotation(finalword, match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing, base, annot_stopper)
                                                                                skipper = True
                                                                                annot_stopper = True
                                                                            # [Any] Only one word, so continue to the next word.   (middle of text)
                                                                            else:
                                                                                self.AddAnnotation(finalword, match, self.count, m, " ", taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing, base, annot_stopper)
                                                                                skipper = False
                                                                                annot_stopper = True
                                                            # [Any] Only one word, so continue to the next word.   (words at the end of the text)
                                                            elif index == len(wordlist) -1 :
                                                                self.AddAnnotation(finalword, match, self.count, m, " ", taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing, base, annot_stopper)
                                                                skipper = False
                                                                annot_stopper = True
                                                    
                                                    else:
                                                            if index <= len(wordlist) -2:
                                                                    finalword = list(finalword)
                                                                    nextword = wordlist[index + 1]
                                                                    nextword = self.RemovePunc(nextword, []) 
                                                                    nonregistered_genus = []
                                                                    possible_species = "".join(finalword) + " " + str(nextword)
                                                                    for cn in CleanNames:
                                                                            cn = cn.split(" ")
                                                                            if "".join(finalword) in cn:
                                                                                    nonregistered_genus.append(" ".join(cn))
                                                                            else: continue  
                                                                    if len(nonregistered_genus) >= 1:
                                                                            for nrg in nonregistered_genus:
                                                                                if nrg == possible_species:
                                                                                        match = next((l for l in dict_data if l['CleanName'] == nrg), None)
                                                                                        modifier = "species"
                                                                                        self.AddAnnotation(possible_species, match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing, base, annot_stopper)
                                                                                        skipper = True
                                                                                        annot_stopper = True  
                                                                    #G. species
                                                                    elif (len(finalword) == 2 and finalword[1] == '.' and finalword[0].isupper) or (len(finalword) == 1 and finalword[0].isupper):
                                                                                for cn in CleanNames:
                                                                                        cn = cn.split(" ")
                                                                                        if nextword in cn and cn[0][0] == finalword[0]:
                                                                                                possible.append(" ".join(cn))
                                                                                        else: continue  
                                                                                
                                                                                if len(possible) == 1:
                                                                                        match = match = next((l for l in dict_data if l['CleanName'] == possible[0]), None)
                                                                                        modifier = "species"
                                                                                        self.AddAnnotation(possible_species, match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing, base, annot_stopper)
                                                                                        skipper = True
                                                                                        annot_stopper = True
                                                                                elif len(possible) >= 2:
                                                                                        generalcheck = True
                                                                                        for tpf in taxa_per_file:
                                                                                                    if tpf in possible:
                                                                                                            match = next((l for l in dict_data if l['CleanName'] == tpf), None)
                                                                                                            possible_species = "".join(finalword) + " " + str(nextword)
                                                                                                            modifier = "species"
                                                                                                            self.AddAnnotation(possible_species, match, self.count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing, base, annot_stopper)
                                                                                                            skipper = True
                                                                                                            annot_stopper = True
                                                                                                            generalcheck == False
                                                                                                    else:
                                                                                                        continue
                                                                                        for p in possible:
                                                                                            match = next((l for l in dict_data if l['CleanName'] == p), None)
                                                                                            duptxids.append(match['TaxID'])
                                                                                        if generalcheck == True:  
                                                                                            strains[possible_species] = len(duptxids)
                                                                                            match = next((l for l in dict_data if l['TaxID'] == duptxids[0]), None) 
                                                                                            self.AddAnnotation(possible_species, match, self.count, m, " ", taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing, base, annot_stopper)
                                                                                            annot_stopper = True
                                                                                            skipper = True
                                                                                        else:
                                                                                            continue
                                                                                else:
                                                                                    continue
                                                                    elif newword in CleanNames: 
                                                        
                                                                        match = next((l for l in dict_data if l['CleanName'] == newword), None)
                                                                        self.AddAnnotation(newword, match, self.count, m, " ", taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing, base, annot_stopper)
                                                                        skipper = False
                                                                        annot_stopper = True 
                                                                    
                                                                    else:
                                                                        continue
                                                            else:
                                                                if newword in CleanNames: 
                                                        
                                                                        match = next((l for l in dict_data if l['CleanName'] == newword), None)
                                                                        self.AddAnnotation(newword, match, self.count, m, " ", taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing, base, annot_stopper)
                                                                        skipper = False
                                                                        annot_stopper = True  
                                                                else:
                                                                    sentenceoffset += (len(word)+ 1)
                                                                    continue
                                                    sentenceoffset += (len(word)+ 1) 
                                                except:
                                                    problemwords.append(finalword, sentenceoffset)
                                        bar()
                                        # Post processing
                                        #Adjust offset
                                        for ann, ian in enumerate(m['annotations']):
                                              wordfound = ian['text']
                                              currentoffset = ian['locations']['offset']
                                              locs = []
                                              distances = []
                                              
                                              for start in re.finditer(wordfound, textsection):
                                                    locs.append({'offset': start.start(), 'distance': abs(currentoffset - (start.start() + m['offset']))})
                                                    #print(wordfound, {'offset': start.start(), 'distance': abs(currentoffset - (start.start() + m['offset']))})
                                              if len(locs) >=2:
                                                    for i in locs: 
                                                        distances.append(i['distance'])
                                                    min_val = min(distances)
                                                    correct = 0
                                                    for i in locs: 
                                                        if i['distance'] == min_val:
                                                            correct =  i['offset']
                                                        ian['locations']['offset'] = int(correct) + int(m['offset'])
                                              elif len(locs) == 1 :
                                                    ian['locations']['offset'] = int(locs[0]['offset']) + int(m['offset'])
                                              else:
                                                    continue
                                        # Adjust for strains
                                        if needs_processing != []:
                                             for ann, ian in enumerate(m['annotations']):
                                                if  '[' in ian['infons']['identifier']: # Finds ambiguous annotations
                                                     offset = ian['locations']['offset']
                                                     countog = ian["id"]
                                                     # Look at surrounding annotations
                                                     others = []
                                                     if ann == 0 and len(m['annotations']) != 1:
                                                        others.append(m['annotations'][ann + 1]['infons']['type'])
                                                     elif ann <= len(m['annotations']) - 2:
                                                        if '[' not in m['annotations'][ann -1]['infons']['identifier']:
                                                            others.append(m['annotations'][ann - 1]['infons']['type'])
                                                        elif '[' not in m['annotations'][ann +1]['infons']['identifier']:
                                                             others.append(m['annotations'][ann + 1]['infons']['type'])
                                                        else: ian['infons']['type'] = 'unresolved'
                                                                
                                                     else:
                                                            ian['infons']['type'] = 'unresolved'
                                                     loc_id = ian['infons']['identifier'].find('[')
                                                     unproc_str = ian['infons']['identifier'][loc_id:len(ian['infons']['identifier'])].split(" ")
                                                     possible_ids = []
                                                     items = []
                                                     for pid in unproc_str: 
                                                          pid = pid.strip("',[]")
                                                          possible_ids.append(pid)
                                                     for pid in possible_ids:
                                                          
                                                          match = next((l for l in dict_data if l['TaxID'] == pid), None)
                                                          identifier = self.MakeIdentifier(match," ", "")  
                                                          item = [pid, identifier]
                                                          items.append(item)
                                                     unresolved  = False
                                                     for item in items:
                                                          if item[1] in others:
                                                               match = next((l for l in dict_data if l['TaxID'] == item[0]), None)
                                                               dictannot = {
                                                                                "text":ian['text'],
                                                                                "infons":{
                                                                                    "identifier": item[0],
                                                                                    "type": item[1] ,
                                                                                    "annotator":"dhylan.patel21@imperial.ac.uk",
                                                                                    "date": time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()) ,
                                                                                    "parent_taxonomic_id": match['ParentTaxID']
                                                                                },
                                                                                "id": countog,
                                                                                "locations":{
                                                                                    "length": len(ian['text']),
                                                                                    "offset": offset ,
                                                                                    
                                                                                }
                                                                            }
                                                               m['annotations'].append(dictannot)            
                                                               unresolved = False
                                                               break
                                                          else:
                                                               unresolved = True
                                                               continue
                                                     if unresolved == True:
                                                          ian['infons']['type'] = 'unresolved'
                                                          
                                        m['annotations'].sort(key = lambda e: (e["id"]))    
                                        for ian in m['annotations']:
                                             if  '[' in ian['infons']['identifier'] and ian['infons']['type'] != 'unresolved':
                                                  m['annotations'].remove(ian)
                                             else:
                                                  continue
                                        m['annotations'].sort(key = lambda e: (int(e["id"])))  
                                        for ian in m['annotations']:
                                              if ian['infons']['type'] == 'unresolved':
                                                 loc_id = ian['infons']['identifier'].find('[')
                                                 unproc_str = ian['infons']['identifier'][loc_id:len(ian['infons']['identifier'])].split(" ")
                                                 possible_ids = []
                                                 itemstoadd = []
                                                 trimmed = 0
                                                 parentids = []
                                                 for pid in unproc_str: 
                                                    pid = pid.strip("',[]")
                                                    possible_ids.append(pid)
                                                 if len(possible_ids) >= 4:
                                                       trimmed = len(possible_ids)
                                                       possible_ids = [possible_ids[0], possible_ids[len(possible_ids) -1]]
                                                 for pid in possible_ids:
                                                    match = next((l for l in dict_data if l['TaxID'] == pid), None)
                                                    identifier = self.MakeIdentifier(match," ","")    
                                                    itemstoadd.append(identifier)
                                                    parentids.append(match['ParentTaxID'])
                                                 ian['infons']['identifier'] = possible_ids
                                                 ian['infons']['type'] = itemstoadd
                                                 ian['infons']['parent_taxonomic_id'] = parentids
                                                 if trimmed != 0:
                                                    ian['infons']['identifier'] = str(possible_ids) + " Trimmed from : " + str(trimmed)
                                                    ian['infons']['type'] = itemstoadd
                                                    ian['infons']['parent_taxonomic_id'] = parentids

                                                 
                                        
                                taxa_per_file = {*taxa_per_file}
                                taxa_per_file = list(taxa_per_file)
                                print(taxa_per_file)
                                print(problemwords)
                                taxalist.append(taxa_per_file)
                                for i in data: 
                                    a_file = open(self.output_directory + folder + "/" +str(in_file), "w")
                                    json.dump(data, a_file, indent = 4)
                                    a_file.close()
             
            stop_time = datetime.datetime.now() #stop time
            message = 'Start time is ' + str(start_time) + '\n' + 'Stop time is ' + str(stop_time)  

            with open(self.output_directory + folder + "/" + 'runtime.log', 'a') as time_file:   #start and stop time are written into a runtime.log file
                
                
                time_file.write(message)

            with open(self.output_directory + folder + "/" + 'taxafound.log', 'a') as taxa_file:   #All taxa found are added to taxafound.log file.
                
                taxa_file.write("This file contains all the kingdoms and taxonomic ranks for each species found.")
                taxa_file.write(str(taxalist))
                taxa_file.write("This file contains all problem words.")
                taxa_file.write(str(problemwords))

        def AddAnnotation(self, word, match, count, m, modifier, taxa_per_file, sentenceoffset, offsetoftext, strains, dict_data, duptxids, idinuse, needs_processing, base, annot_stopper):
            if match['TaxID'] not in idinuse and annot_stopper == False:
                 
                self.count = int(count) + 1
                if word in strains:
                    
                    #repeats = int(strains[word])
                    if duptxids == []:
                        for i in dict_data:
                            if i['CleanName'] == word:
                                duptxids.append(i['TaxID'])
                   
                    typels = []
                    parentls = []
                      
                    for i in duptxids:
                            
                                match = next((l for l in dict_data if l['TaxID'] == i), None)
                                
                                typemod = self.MakeIdentifier(match, " ", "") 
                                typels.append(typemod)
                                parentls.append(match['ParentTaxID'])
                if  duptxids == []:
                     identifierstring = match['TaxID']
                elif modifier != " ":
                     identifierstring = match['TaxID']
                else: 
                     identifierstring = str(duptxids)
                     needs_processing.append(identifierstring)
                if modifier != " ": 
                        dictannot = {
                                            "text":word,
                                            "infons":{
                                                "identifier": identifierstring,
                                                "type": self.MakeIdentifier(match, modifier, "") ,
                                                "annotator":"dhylan.patel21@imperial.ac.uk",
                                                "date": time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()) ,
                                                "parent_taxonomic_id": match['ParentTaxID']
                                            },
                                            "id": str(base) + str(count),
                                            "locations":{
                                                "length": len(word),
                                                "offset": sentenceoffset + offsetoftext ,
                                                
                                            }
                                        }
                        taxa_per_file.append(str(match['CleanName']))
                elif '[' in identifierstring:
                        
                        dictannot = {
                                            "text":word,
                                            "infons":{
                                                "identifier": identifierstring,
                                                "type": typels ,
                                                "annotator":"dhylan.patel21@imperial.ac.uk",
                                                "date": time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()) ,
                                                "parent_taxonomic_id": parentls
                                            },
                                            "id": str(base) + str(count),
                                            "locations":{
                                                "length": len(word),
                                                "offset": sentenceoffset + offsetoftext ,
                                                
                                            }
                                        }
                        taxa_per_file.append(word)
                else:
                        dictannot = {
                                            "text":word,
                                            "infons":{
                                                "identifier": identifierstring,
                                                "type": self.MakeIdentifier(match, modifier, "")  ,
                                                "annotator":"dhylan.patel21@imperial.ac.uk",
                                                "date": time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()) ,
                                                "parent_taxonomic_id": match['ParentTaxID']
                                            },
                                            "id": str(base) + str(count),
                                            "locations":{
                                                "length": len(word),
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
        def RemovePunc (self, word, finalword):
                wordaslist = []
                listword = list(word)
                for i in listword:
                      if i == '(' or i == ')':
                            continue
                      else:
                          wordaslist.append(i)
                if len(wordaslist)!=2:   
                    for i in wordaslist:
                        if i == '.' or i == ',' or i == '(' or i == ')':
                            continue
                        else:
                            finalword.append(i)
                    
                elif wordaslist[0].isupper() and wordaslist[1] == '.':
                         for i in wordaslist:
                                    finalword.append(i)
                else:
                     for i in wordaslist:
                        if i == '.' or i == ',' or i == '(' or i == ')':
                            continue
                        else:
                            finalword.append(i)
                     
                finalword = "".join(finalword)
                return finalword
        def MakeIdentifier (self, match, modifier, identifier):
              kingdic = {'NCBI:txid2': 'bacteria', 'NCBI:txid2157':'archaea', 'NCBI:txid4751':'fungi'}
              kingdom = kingdic[match['KingdomID']]
              if modifier != " ":
                    identifier = kingdom + "_" + modifier
              else:
                    identifier = kingdom + "_" + match['TaxRank']
              return identifier
        

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
