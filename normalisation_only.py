import json
from collections import Counter
# Most up to date version of this dictionary. 
opendic = open('NCBI_tax_dictionary10.json', encoding = 'utf-8') 
dict_data = json.load(opendic)
CleanNames = [i['CleanName'] for i in dict_data]
#Used to produce the strains dictionary for duplicate taxa where clean names have multiple mentions. 
cncounts = dict(Counter(CleanNames))
strains = {key:value for key, value in cncounts.items() if value > 1}
CleanNames =list(dict.fromkeys(CleanNames))
lcCleanNames = [i['CleanName'].lower() for i in dict_data]
# Used to demonstrate that the program is working. They can be removed in place of actual data.
testwords = ["Escherichia coli", "Clostridium botulinum", "Inventus fakus"]
for i in testwords:
    if i in CleanNames: # Adjust to lcCleanNames if only using lower case values

        match = next((l for l in dict_data if l['CleanName'] == i), None)
        normalised_identifier = match["TaxID"]
        print(normalised_identifier)
    else:
        print("Name not recognised")