import requests
import time
import os
class Childnodes:
        def __init__(self, childnodesdir, specfound, ncbikey):
            #Initialise inputs
            self.childnodesdir = childnodesdir
            self.specfound = specfound
            self.ncbikey = ncbikey
        def updatenewspec (self):
                with open(self.childnodesdir, "a+", encoding = 'utf-8') as fp:
                                    text2 = fp.readlines()
                                    firsttry = []
                                    
                                    for i in text2:
                                          id = i.split(" ")[0]
                                          firsttry.append(str(id))
                                    for sf in self.specfound:
                                            if sf in firsttry:
                                                    continue
                                            else:
                                                    url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/" + str(sf) + "/filtered_subtree"
                                                    response = requests.request("GET", url, params =  {"key": self.ncbikey})

                                                    
                                                    if response.status_code == 200:
                                                            jsonresponse = response.json()
                                                            if len(jsonresponse["edges"].keys()) >=2:
                                                                actualkeys = jsonresponse["edges"].keys()
                                                                fp.write(sf+ " " + str(len(actualkeys) - 2 ) + "\n")
                                                    else:  
                                                            fp.write(sf+ " " + "\n")
        def filerefresh (self):
                with open(self.childnodesdir, "r", encoding = 'utf-8') as fp:
                                    text2 = fp.readlines()
                                    firsttry = []
                                    
                                    for i in text2:
                                          id = i.split(" ")[0]
                                          firsttry.append(str(id))
                filename = str(time.strftime('%Y-%m-%d_%H-%M-%S')) +'childnodes.txt'
                with open(filename, "a", encoding = 'utf-8') as fp:
                        for ft in firsttry:
                                url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/" + str(ft) + "/filtered_subtree"
                                response = requests.request("GET", url , params =  {"key": self.ncbikey})
                                if response.status_code == 200:
                                        jsonresponse = response.json()
                                        if len(jsonresponse["edges"].keys()) >=2:
                                            actualkeys = jsonresponse["edges"].keys()
                                            fp.write(ft+ " " + str(len(actualkeys) - 2 ) + "\n")
                                else:
                                        
                                        fp.write(ft+ " " + "\n")
                        fp.close()
                with open(str(time.strftime('%Y-%m-%d_%H-%M-%S')) +'childnodes.txt', "r", encoding = 'utf-8') as fp:
                        
                        linecount =  len(fp.readlines())
                if os.path.isfile(filename) == True and linecount == len(firsttry):

                    os.remove(self.childnodesdir)