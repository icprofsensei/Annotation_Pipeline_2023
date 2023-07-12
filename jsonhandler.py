#json handler
import json
import csv
with open ("common words.csv", newline = '') as g:
    reader = csv.reader(g)
    datalist = list(reader)
    print(datalist)
    
with open("NCBI_TAX_DIC.json", 'r+') as f:
    data = json.load(f)
    newdata = []
    cleannames = []
    for i in data:    
        if i["CleanName"] in datalist[0]:
            print("found_one")
        elif i["CleanName"] not in cleannames:

            newdata.append(i)
            cleannames.append(i["CleanName"])
with open("newdic.json", "w") as e:
    json.dump(newdata, e, indent = 4)