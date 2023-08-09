#tsv_cleaner_to_json
import json
import csv
from alive_progress import alive_bar
import os
array = []
tsv = open('NCBI_taxonomy_dictionary_v8.tsv', 'r')
titles = ["OriginalName", "TaxRank", "CleanName", "TaxID", "ParentTaxID", "KingdomID"]
for line in tsv:
    d = {}
    #split each line when the tab is found. t relates to titles and f relates to the line being split. They are added to d. 
    for t, f in zip(titles, line.split('\t')):
        d[t] = f.strip()
    array.append(d)
with open("NCBI_TAX_DIC.json", 'w', encoding='utf-8') as output_file:
    output_file.write(json.dumps(array, indent = 4))

#json handler

with open ("common_words.csv", newline = '') as g:
    reader = csv.reader(g)
    datalist = list(reader)
    print(datalist)
    
with open("NCBI_TAX_DIC.json", 'r+') as f:
    data = json.load(f)
    newdata = []
    with alive_bar(len(data)) as bar:
        for i in data:    
            if i["CleanName"] in datalist[0]:
                print(i)
            else:
                newdata.append(i)
            bar()

with open("NCBI_tax_dictionary8.json", "w") as e:
    json.dump(newdata, e, indent = 4)

if os.path.isfile('NCBI_tax_dictionary8.json') == True:
                 os.remove('NCBI_TAX_DIC.json')