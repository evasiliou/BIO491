import pandas as pd
import csv
from UniProtMapper 

def uniprot_id_mapping (query_list, uniprot_acc): #in case that the uniprot id mapper does not work
  '''
    query_list: list of uniprot ids
    uniprot_acc: bool value
  '''
  if uniprot_acc == True:
    filename = 'idmapping_acc.csv'
  else:
    filename = 'idmapping_id.csv'

  df = pd.read_csv (filename)
  uniprot_ids = list (df['Uniprot ID'])
  gene_names = list (df['Gene Name'])

  # print (uniprot_ids)
  # print (gene_names)

  not_found = []
  dictionary = {}

  # print (query_list)

  for protein in query_list:
    if protein in uniprot_ids:
      index = uniprot_ids.index(protein)
      # print (index)
      gene = gene_names[index]
      dictionary[protein] = {'Gene': gene}

    else:
      dictionary[protein] = {'Gene': 'None'}
      not_found.append(protein)
      with open ('not_found_proteins.txt', 'a') as file: #append the proteins that were not found in the 'not_found_proteins.tsv' file
        if protein not in 'not_found_proteins.txt':
          file.write(f'%s\n' % (protein))


  print ('\n Converted Proteins:', len(dictionary))
  print (dictionary)

  print ('\n Not found:', len(not_found))


  return dictionary #returns a dictionary
