
import os
import json
import requests
class species_miner:
        def __init__(self, folder_directory, keyword):
            #Initialise inputs
            self.folder_directory = folder_directory
            self.keyword = keyword
        
        def Miner(self):

                all_files = os.listdir(self.folder_directory)   
                PMC_files = [n for n in all_files if n.endswith('bioc.json')]
                annotation_list = []
                with open('txidlist.txt', 'w') as f: # self.folder_directory+ "/" + 
                    for i in PMC_files:
                        in_file = i
                        #print("Annotating file: ",  in_file) 
                        with open(self.folder_directory+ "/" + in_file , encoding = 'utf-8') as m_file:
                            data = json.load(m_file)
                            documents=data['documents'] 
                            for j in documents:
                                    passages= j['passages'] 
                                    for m in passages:
                                        self.keyword = str(self.keyword)
                                        important = False
                                        upperword = self.keyword[0].upper() + self.keyword[1:len(self.keyword)]
                                        lowerword = self.keyword[0].lower() + self.keyword[1:len(self.keyword)]
                                        for v in m['infons'].values():
                                                        if v == upperword or v == lowerword:
                                                            important = True
                                                        elif upperword in v.split(" "):
                                                            important = True
                                                        elif lowerword in v.split(" "):
                                                            important = True
                                                        else:
                                                            continue
                                        
                                    
                                        if important == True or self.keyword == 'ALL':
                                                
                                            for ian in m['annotations']:
                                                annotation_list.append(str(ian['infons']['identifier']))
                                        else:
                                            continue
                        annotation_list = list(dict.fromkeys(annotation_list))
                        annotation_list = [i.lstrip("NCBI:txid") for i in annotation_list if i.startswith('[') == False]
                        querystring = ""
                        for i in annotation_list:
                              querystring += i
                              querystring += "%2C"
                        url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/" + querystring
                        response = requests.request("GET", url)
                        if response.status_code == 200:
                                    print("sucessfully fetched the data")
                                    
                        else:
                                    print(f"There's a {response.status_code} error with your request")
                        jsonresponse = response.json()
                        for i in jsonresponse['taxonomy_nodes']:
                            #print(i)
                            if i['query'] == ['']:
                              continue
                            else:
                              f.write("Query")
                              f.write(str(i['query']))
                              f.write("Lineage")
                              f.write(str(i['taxonomy']['lineage']))
                              f.write('\n')
                    #print(annotation_list)                      

        def Analyser(self):
                          with open('txidlist.txt',encoding = 'utf-8') as m_file:
                                text = m_file.readlines()
                                size = len(text)
                                print(size)
                                empty = []
                                for line in text:
                                      if line not in empty:
                                            empty.append(line)
                                output_file = "output.txt"
                                with open(output_file, "w") as fp:
                                    fp.write("".join(empty))
                          with open(self.folder_directory + 'outputclean.txt',encoding = 'utf-8') as fp:
                                    text2 = fp.readlines()
                                    size2 = len(text2)
                                    print(size2)
                          all_files = os.listdir(self.folder_directory)   
                          PMC_files = [n for n in all_files if n.endswith('bioc.json')]
                          annotation_list = []
                          for i in PMC_files:
                                in_file = i
                                with open(self.folder_directory+ "/" + in_file , encoding = 'utf-8') as m_file:
                                    data = json.load(m_file)
                                    documents=data['documents'] 
                                    for j in documents:
                                            passages= j['passages'] 
                                            for m in passages:
                                                self.keyword = str(self.keyword)
                                                important = False
                                                upperword = self.keyword[0].upper() + self.keyword[1:len(self.keyword)]
                                                lowerword = self.keyword[0].lower() + self.keyword[1:len(self.keyword)]
                                                for v in m['infons'].values():
                                                                if v == upperword or v == lowerword:
                                                                    important = True
                                                                elif upperword in v.split(" "):
                                                                    important = True
                                                                elif lowerword in v.split(" "):
                                                                    important = True
                                                                else:
                                                                    continue
                                                
                                            
                                                if important == True or self.keyword == 'ALL':
                                                        
                                                    for ian in m['annotations']:
                                                        annotation_list.append(str(ian['infons']['identifier']))
                                                else:
                                                    continue
                          annotation_list = list(dict.fromkeys(annotation_list))
                          annotation_list = [i.lstrip("NCBI:txid") for i in annotation_list if i.startswith('[') == False]
                          print(len(annotation_list))
                          with open(self.folder_directory +'outputclean.txt',encoding = 'utf-8') as fp:
                                    text2 = fp.readlines()
                                    firsttry = []
                                    for i in text2:
                                          id = i[6:i.find(']')]
                                          id = id.replace("'", "")
                                          firsttry.append(str(id))
                          print(firsttry)
                          secondtry = []
                          for id in annotation_list:
                                if id in firsttry:
                                      continue
                                else:
                                      secondtry.append(id)
                          print(secondtry)
                          with open('txidlist.txt', 'a') as f:
                                querystring = ""
                                for i in secondtry:
                                    querystring += i
                                    querystring += "%2C"
                                url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/" + querystring
                                response = requests.request("GET", url)
                                if response.status_code == 200:
                                            print("sucessfully fetched the data")
                                            
                                else:
                                            print(f"There's a {response.status_code} error with your request")
                                jsonresponse = response.json()
                                for i in jsonresponse['taxonomy_nodes']:
                                    #print(i)
                                    if i['query'] == ['']:
                                        continue
                                    else:
                                        f.write("Query")
                                        f.write(str(i['query']))
                                        f.write("Lineage")
                                        f.write(str(i['taxonomy']['lineage']))
                                        f.write('\n')