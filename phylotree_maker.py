import json

class TreeMaker:
        def __init__(self, annotated_file_path, NCBI_tax_dic):
            #Initialise inputs
            self.annotated_file_path = annotated_file_path
            self.NCBI_tax_dic = NCBI_tax_dic

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
            
        def dictionary_to_lineage(self, text):
             opendic = open(self.NCBI_tax_dic, encoding = 'utf-8') 
             dict_data = json.load(opendic)
             Tax_parent = [{'TaxID': i['TaxID'], 'Rank': i['TaxRank'], 'ParentTaxID': i['ParentTaxID']} for i in dict_data]
             
             lineage = []
             lineage.append(text)
             
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
                elif lineagemaker == False:
                        print(lineage)
                        
                     
             
                                                