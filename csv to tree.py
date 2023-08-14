import json
import csv
from alive_progress import alive_bar
opendic = open('NCBI_tax_dictionary8.json', encoding = 'utf-8') 
dict_data = json.load(opendic)
species = [i['TaxID'] for i in dict_data if i['TaxRank'] == 'species']
parentids = [i['ParentTaxID'] for i in dict_data if i['TaxRank'] == 'species']
species = list(dict.fromkeys(species))
parentids = list(dict.fromkeys(parentids))
with open('Species_lineages.csv', mode='w') as species_file:
    species_writer = csv.writer(species_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    with alive_bar(len(species)) as bar:
        for i in species:
                    speclist = self.dictionary_to_lineage(i)
                    if len(speclist) == 0:
                        continue
                    else:
                        species_writer.writerow(speclist)
                    bar()