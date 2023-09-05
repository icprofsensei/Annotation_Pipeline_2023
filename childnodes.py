import requests
with open('output.txt',encoding = 'utf-8') as fp:
                                    text2 = fp.readlines()
                                    firsttry = []
                                    output_file = "childnodes.txt"
                                    
                                    for i in text2:
                                          id = i[6:i.find(']')]
                                          id = id.replace("'", "")
                                          firsttry.append(str(id))
                                    for i in firsttry:
                                        url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/" + str(i) + "/filtered_subtree"
                                        response = requests.request("GET", url)
                                        if response.status_code == 200:
                                                #print("sucessfully fetched the data")
                                                jsonresponse = response.json()
                                                #print(jsonresponse["edges"][i])
                                                if len(jsonresponse["edges"].keys()) >=2:
                                                    keys = jsonresponse["edges"].items()
                                                    actualkeys = jsonresponse["edges"].keys()
                                                    with open(output_file, "a") as fp:
                                                                    fp.write(i+ " " + str(len(actualkeys) - 2 ) + "\n")
                                                
                                                    '''
                                                    if i in actualkeys:
                                                            with open(output_file, "a") as fp:
                                                                    fp.write(i+ " " + str(len(jsonresponse["edges"][i]['visible_children'])) + "\n")
                                                    else:
                                                                
                                                        second = list(keys)[1]
                                                        #print(list(second)[1])
                                                        with open(output_file, "a") as fp:
                                                            fp.write(i+ " " + str(len(second[1]['visible_children'])) + "\n")
                                                            '''
                                        else:
                                                print(f"There's a {response.status_code} error with your request")
                                                with open(output_file, "a") as fp:
                                                            fp.write(i+ " " + "\n")
                                        