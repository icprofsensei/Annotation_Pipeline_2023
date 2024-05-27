import json
import re
 
def CheckLatin (word, newword):
        if word.endswith('a') and not word.endswith('acteria'):
            newword = re.sub('a$', 'ae', word)
        elif word.endswith('us'):
            newword = re.sub('us$', 'i', word)
        elif word.endswith('um'):
            newword= re.sub('um$', 'a', word)
        elif word.endswith('is'):
            newword= re.sub('is$', 'es', word)  
        elif word.endswith('e'):
            newword= re.sub('e$', 'es', word)    
        else:
             newword = word
        return newword
# Opening JSON file
f = open('New_Pipeline_2023/NCBI_tax_dictionary11_abb.json')
 
# returns JSON object as 
# a dictionary
data = json.load(f)
 
# Iterating through the json
# list
newdata = []
for i in data:
    OG = i['OriginalName']
    latpl = CheckLatin(str(OG), '')
    #print(OG, latpl)
    i['Latin_plural'] = latpl
    newdata.append(i)
# Closing file
f.close()
json_string  = json.dumps(data)
with open("New_Pipeline_2023/NCBI_tax_dictionary11_abb_latpl.json", "w") as e:
    json.dump(data, e, indent = 4)