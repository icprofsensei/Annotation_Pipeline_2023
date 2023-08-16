import json

class TreeMaker:
        def __init__(self, NCBI_tax_dic, text):
            #Initialise inputs
            #self.annotated_file_path = annotated_file_path
            self.NCBI_tax_dic = NCBI_tax_dic
            self.text = text

        def dictionary_to_lineage(self):
             opendic = open(self.NCBI_tax_dic, encoding = 'utf-8') 
             dict_data = json.load(opendic)
             Tax_parent = [{'TaxID': i['TaxID'],'ParentTaxID': i['ParentTaxID']} for i in dict_data]
             
             lineage = []
             lineage.append(self.text)
             
             lineagemaker = True
             for i in lineage:
                if lineagemaker == True:
                
                    for tp in Tax_parent:
                            if tp['TaxID'] == i:
                                lineage.append(tp['ParentTaxID'])
                                
                                break
                    
                    if i in ['NCBI:txid4751', 'NCBI:txid2157', 'NCBI:txid2']:
                        lineagemaker = False
                    else:
                        lineagemaker = True
                else:
                        return(lineage)
             return(lineage)
                        
  
                           
             
'''                              
        def TM(self):
            with open(self.annotated_file_path , encoding = 'utf-8') as m_file:
                taxids_in_article = []
                data = json.load(m_file)
                documents=data['documents']
                for j in documents:
                            passages= j['passages'] 
                            for m in passages: 
                                    annotations = m['annotations'] 
                                    if len(annotations) > 0:
                                            for ann in annotations:
                                                    taxids_in_article.append(ann['infons']['identifier'])
            for i in taxids_in_article:
                  self.dictionary_to_lineage(i)
        

        def Tree(self):
              opendic = open(self.NCBI_tax_dic, encoding = 'utf-8') 
              dict_data = json.load(opendic)
              species = [i['TaxID'] for i in dict_data if i['TaxRank'] == 'species']
              species = list(dict.fromkeys(species))
              print(species)
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
'''  