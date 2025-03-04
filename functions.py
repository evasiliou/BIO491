import pandas as pd
import csv
import numpy as np
from gprofiler import GProfiler
from ast import literal_eval
import requests
import igraph as ig
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles
from UniProtMapper import ProtMapper
import time
from bs4 import BeautifulSoup
import re
import sklearn as sk
from sklearn import linear_model
from sklearn.linear_model import LinearRegression
import scipy.stats as stats
from scipy.stats import mannwhitneyu
from math import log
from math import log10


not_found_file = 'not_found_proteins.txt' #create a file that will save the proteins not found
with open (not_found_file, 'w') as file: #just to create the file
  file.write('')

def fetching (total_size):

  # Define total number of items and batch size
  total_items = total_size
  batch_size = 500
  fetched_items = 0  # Counter to track progress

  while fetched_items < total_items: #i tried to run while fetched_items <= total_items, in order to take the last item of the total proteins (e.g. all_proteins[:total_items])
  #but i fall into an infinite loop
      # Calculate the next batch size
      remaining_items = total_items - fetched_items
      current_batch_size = min(batch_size, remaining_items)

      # Simulate the data fetching (replace with your fetch function)
      # print(f"Fetching batch: {fetched_items + 1} to {fetched_items + current_batch_size}...")
      time.sleep(0.5)  # Simulating delay for each fetch; remove this line in actual fetching

      # Update the count of fetched items
      fetched_items += current_batch_size

      # Calculate and display progress
      progress_percentage = (fetched_items / total_items) * 100
      # print(f"Progress: {fetched_items}/{total_items} ({progress_percentage:.2f}%)\n")

  print("\n All items fetched!")


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

def protein_to_protein_accession_conversion (query_list):

    mapper = ProtMapper()

    query = ','.join(query_list)

    from_db = 'UniProtKB_AC-ID'
    to_db = 'UniProtKB'
    ids = query


    results_df, failed_list  = mapper.get(ids=ids, from_db=from_db, to_db=to_db) #mapper.get returns a tuple with a df of the results, and a list with the ids that were not found
    #just a note: the failed list contains all the ids in a single string e.g ['P49913,P27695,P0DMM9']

    # fetching(len(query_list))

    dictionary = {}
    for i in range (len(results_df)):
      dictionary[results_df.loc[i, 'From']] = {'Acc': results_df.loc[i, 'Entry']}

    print ('\n Converted Proteins:', len(dictionary))
    print (dictionary)

    element = failed_list[0] #the element is a string that contains all the query proteins (e.g. 'P49913,P27695,P0DMM9')
    all_proteins = element.split(',') #this is a list


    # print (all_proteins)
    print ('\n All proteins:', len(all_proteins))

    common_proteins = list (set(all_proteins).intersection(set(dictionary.keys())))
    print ('\n Converted Proteins:', len(common_proteins))

    not_found= []

    for protein in all_proteins:
      if protein not in common_proteins:
        not_found.append(protein)
        dictionary[protein] = {'Acc': 'None'}
        with open ('not_found_proteins_acc.txt', 'w') as file: #append the proteins that were not found in the 'not_found_proteins.tsv' file
          if protein not in 'not_found_proteins_acc.txt':
            file.write(f'%s\n' % (protein))

    print ('\n Not found:', len(not_found))


    return dictionary #returns a dictionary

def protein_to_gene_conversion(query_list):
    '''
    The protein_to_gene_conversion function takes as argument a list with uniprot ids and returns a dictionary
    that has as key the uniprot id and as value a dictionary with equivalent gene name {'Uniprot id': {'Gene': Gene name}}.
    '''
    mapper = ProtMapper()

    query = ','.join(query_list)

    from_db = 'UniProtKB_AC-ID'
    to_db = 'Gene_Name'
    ids = query


    results_df, failed_list  = mapper.get(ids=ids, from_db=from_db, to_db=to_db) #mapper.get returns a tuple with a df of the results, and a list with the ids that were not found
    #just a note: the failed list contains all the ids in a single string e.g ['P49913,P27695,P0DMM9']

    # fetching(len(query_list))

    dictionary = {}
    for i in range (len(results_df)):
      dictionary[results_df.loc[i, 'From']] = {'Gene': results_df.loc[i, 'To']}

    print ('\n Converted Proteins:', len(dictionary))
    print (dictionary)

    element = failed_list[0] #the element is a string that contains all the query proteins (e.g. 'P49913,P27695,P0DMM9')
    all_proteins = element.split(',') #this is a list


    # print (all_proteins)
    print ('\n All proteins:', len(all_proteins))

    common_proteins = list (set(all_proteins).intersection(set(dictionary.keys())))
    print ('\n Converted Proteins:', len(common_proteins))

    not_found= []

    for protein in all_proteins:
      if protein not in common_proteins:
        not_found.append(protein)
        dictionary[protein] = {'Gene': 'None'}
        with open ('not_found_proteins.txt', 'a') as file: #append the proteins that were not found in the 'not_found_proteins.tsv' file
          if protein not in 'not_found_proteins.txt':
            file.write(f'%s\n' % (protein))

    print ('\n Not found:', len(not_found))


    return dictionary #returns a dictionary

def gene_to_protein_conversion(query_list):

    '''
    The gene_to_protein_conversion function takes as argument a list with gene names and returns a dictionary
    that has as key the gene name and as value a dictionary with equivalent protein id {'Gene name': {'Protein': Uniprot id}}.
    '''
    count_none = 0
    gp = GProfiler(
        user_agent='ExampleTool',  # optional user agent
        return_dataframe=True  # return pandas dataframe or plain python structures
    )
    print("\n Total proteins before conversion: ", len(query_list))

    # convert the uniprot ids to the equivalent gene accession numbers. the results are saved in a df
    # query_genes = ','.join(query_list)

    converted_df = gp.convert(organism='hsapiens',
                              query=query_list,
                              target_namespace='UNIPROTSWISSPROT_ACC')

    print(converted_df.head().to_string())
    # print(converted_df)

    # fetching (len(query_list))
    dictionary = {}

    for i in range(len(converted_df)):
        gene = converted_df.loc[i, 'incoming']
        protein = converted_df.loc[i, 'converted']

        if protein != 'None':
            if gene not in dictionary.keys():
                dictionary[gene] = {'Protein': protein}
        else:
            dictionary[gene] = {'Protein': 'None'}
            count_none += 1
            with open ('not_found_proteins.txt', 'a') as file: #append the proteins that were not found in the 'not_found_proteins.tsv' file
              if protein not in 'not_found_proteins.txt':
                file.write(f'%s\n' % (protein))


    incoming_proteins = list(converted_df.loc[:, 'incoming'])
    print ('\n Incoming_proteins:', len(incoming_proteins))
    print("\n Total proteins after conversion: ", (len(dictionary.keys())-count_none))
    print ("\n Total unconverted proteins: ", count_none)
    print ('Dictionary', dictionary)


    return dictionary #returns a dictionary

def gene_enrichment(dictionary, filename):

    '''
    The gene_enrichment function takes as arguments a dictionary (describe what this dictionary contains) and the name of the file where
    the results of the gene enrichment analysis will be saved.
    '''

    #query_genes = protein_to_gene_conversion(dictionary) #this is a list
    query_genes = []
    for uniprot_id in dictionary.keys():
        gene = dictionary[uniprot_id]['Gene']
        query_genes.append(gene)

    gp = GProfiler(
        user_agent='ExampleTool',  # optional user agent
        return_dataframe=True  # return pandas dataframe or plain python structures
    )

    query = gp.profile(organism='hsapiens',
                       query=query_genes,
                      #  sources=['GO:MF'],
                       sources=["GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "WP", "HPA", "CORUM", "HP"],
                       # sources=["GO:MF", "GO:CC", "GO:BP"],
                       no_evidences=False)
    #
    print("\n Total genes after enrichment: ", len(query))
    #print("\n Enrichment Results: \n", query.head().to_string())
    # # print (dictionary)
    # #print(list(query.columns.values))
    #
    go_terms = list((query.loc[:, "native"]))  # pathways with their unique accession number
    print('\n Total GO terms:', len(go_terms))
    name = list((query.loc[:, 'name']))  # a description of the pathway
    print('\n Total names:', len(name))
    p_value = list(query.loc[:, "p_value"])
    print('\n Total p-values:', len(p_value))
    intersections = list(query.loc[:, "intersections"])  # a list with the group of proteins that participate in the equivalent pathway
    print('\n Total intersections:', len(intersections))
    #print(intersections)  # prints a list of lists
    #
    # print("\n Total proteins after enrichment: ", len(dictionary))
    # # print (dictionary)
    column_labels = ['p-value', 'GO terms', 'Pathway(name)', 'Intersections'] #this is a list that contains column names for the csv file
    column_data = [p_value, go_terms, name, intersections] #this is a list that contains all the lists with the data that I want to include in the csv file
    new_csv_file(column_labels, column_data, filename, True)

def read_gene_enrichment_results(filename):
    '''
    The read_gene_enrichmenet_results function takes as argument the name of the file that contains the results from
    the gene enrichment analysis and returns a dictionary: 'go term' = {'p-value': p_value, 'Pathway(name)': pathway, 'Intersections': intersections}
    '''

    csv_file = filename
    df = pd.read_csv(csv_file, converters={'Intersections': literal_eval})  # the converter here is applied for the data under the 'Intersections' column
    # since they were saved as str and not as lists
    print(df.head().to_string())
    #print(df.columns.values)
    go_terms = list(df.loc[:, 'GO terms'])

    dictionary = {}
    for i in range(len(go_terms)):
        go_term = go_terms[i]
        p_value = float(df.loc[i, 'p-value'])
        pathway = (df.loc[i, 'Pathway(name)'])
        # print(cellular_process)
        # autophagy = bool(df.loc[i, 'Autophagy'])
        intersections = (df.loc[i, 'Intersections'])
        # print(intersection)
        dictionary[go_term] = {'p-value': p_value, 'Pathway(name)': pathway, 'Intersections': intersections}

    return (dictionary)
    # return (go_terms, p_value, cellular_process, intersections)

def new_csv_file (column_labels, column_data, filename, index_status):

    '''
    The new_csv_file function aims to create a new csv file using a dictionary which is then converted to a df.
    It takes as arguments, a list with the names of the column in the new file, a list of lists with the data under each column, and the name of the file.
    The index_status takes True or False boolean values and the user selects if the column index will be included or not at the new file.
    '''
    size = len(column_data[0])
    dictionary = {}
    for label, data in zip(column_labels, column_data):
        #print (data)
        dictionary[label] = data

    df=pd.DataFrame(dictionary, index = range(size))
    # print (df)

    df.to_csv(filename, index=index_status, sep=',')

    print (f'\n The file {filename} has been saved!')

#the hippie_query function can be a backup in case that the hippie_api does not work
# def hippie_query (query_list,confidence_score):

#     '''
#     This function takes as argument a file that contains only the enriched genes (previously created). From this file takes the uniprot ids and after making the query
#     in hippie it finds the interactions of these proteins that were found experimentally. Then it prints a URL which when is clicked it downloads a tsv file with the results.
#     It also saves the results from hippie in a tsv file 'hippie_results_0.63.tsv' where '0.63' is the confidence score that was used for the analysis.
#     '''


#     hippie_interactions_df = pd.read_table('HIPPIE_all_db.tsv' ) #saves in a df all the interactions in hippie
#     # print (hippie_interactions_df)

#     hippie_query_pairs = [] #this is a list with tuples that has the pairs of all the hippie interactions
#     interactor_A = []
#     interactor_B = []
#     confidence_values_list = []

#     for i in range (len(hippie_interactions_df)):
#       confidence_value = float(hippie_interactions_df.loc[i, 'Confidence Value'])
#       # print (type(confidence_value))
#       interactor_a = hippie_interactions_df.loc[i, 'Gene Name Interactor A']
#       interactor_b = hippie_interactions_df.loc[i, 'Gene Name Interactor B']

#       if confidence_value >= float(confidence_score) and (interactor_a in query_list or interactor_b in query_list):
#         pair = (interactor_a, interactor_b)
#         if pair not in hippie_query_pairs:
#           hippie_query_pairs.append(pair)

#       interactor_A.append(interactor_a)
#       interactor_B.append(interactor_b)
#       confidence_values_list.append(confidence_value)

#     print ('\n Total HIPPIE interactions: ', len(hippie_query_pairs))
#     # print (len(hippie_query_pairs))


#     column_labels = ['Interactor A', 'Interactor B', 'Confidence Value']
#     column_data = [interactor_A, interactor_B, confidence_values_list]
#     filename = 'hippie_results.csv'
#     new_csv_file(column_labels, column_data, filename, True)


#     return hippie_query_pairs #returns a list of tuples that contains the results from the query
def uniprot_query (entries_list, fields_list):

  '''
    entries_list: a list with the protein id (e.g. 'WIPF1_HUMAN'), protein accession (e.g. 'O43516' ), or gene name (e.g. 'WIPF1') of the proteins of interest

    fields_list: is a list with the fields of each entry you want to query in string format (e.g. 'protein_name', 'cc_function'). The field names can be found in the uniprot documentation (https://www.uniprot.org/help/return_fields)

    Just a note regarding the field_list: The following code block runs only with filds_lists 'protein_name' and 'cc_function'. Thus, if any other fields are required, the appropriate changes must be made.
  '''
  query_list = ','.join(entries_list)

  url = 'https://rest.uniprot.org/uniprotkb/search?'
  
  error_messages = {
      '400':' Bad Request: The request was invalid.',
      '401': 'Unauthorized: Authentication failed.',
      '404': 'Not Found: The requested resource does not exist.',
      '414': 'Request-URI Too Long: The request URI is too long.',
      '500': 'Internal Server Error: The server encountered an error.'
  }

  entries_dictionary = {}
  for entry in entries_list:
    params = {
      'query': entry,
      'fields': fields_list
      }
  
    response = requests.get(url, params=params)
    query = response.json()
    if response.status_code == 200:
  
      result = query['results'][0] # query['results'] is a list of dictionaries, so by using the index  i can access the content of the list 
      # print (result)
      protein_name = result['proteinDescription']['recommendedName']['fullName']['value']
      print (f'\n Protein name: {protein_name}')
      protein_description = result['comments'][0]['texts'][0]['value'] #comments and texts are lists with dictionaries, so by taking the first index i take all the dictionaries out of the list
      print (f'Protein description: {protein_description}')

      entries_dictionary.update({entry:{'Protein Name': protein_name, 'Protein Description': protein_description}})

    elif response.status_code in error_messages.keys():
      print ('\n %s: %s'% (response.status_code, error_messages[str(response.status_code)]))
  
    else:
      print ('\n Error:', response.status_code)
  
  return entries_dictionary 

def hippie_query_api (query_list,confidence_score, filename):

    '''
    This function takes as argument a file that contains only the enriched genes (previously created). From this file takes the uniprot ids and after making the query
    in hippie it finds the interactions of these proteins that were found experimentally. Then it prints a URL which when is clicked it downloads a tsv file with the results.
    It also saves the results from hippie in a tsv file 'hippie_results_0.63.tsv' where '0.63' is the confidence score that was used for the analysis.
    '''

    #Just a note: the beautiful soup package was used with the help of ChatGPT
    url = 'http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/queryHIPPIE.php?'

    error_messages = {
        '400':' Bad Request: The request was invalid.',
        '414': 'Request-URI Too Long: The request URI is too long.',
        '401': 'Unauthorized: Authentication failed.',
        '404': 'Not Found: The requested resource does not exist.',
        '500': 'Internal Server Error: The server encountered an error.'
    }
    # print (query_list)

    query=';'.join(query_list)

    parameters = {
        'proteins':query,
        'layers': '1', #layers = 0 to query interactions within the input set or 1 to query interactions between the input set and HIPPIE (optional, default = 1)
        'conf_thres': confidence_score,
        'out_type': 'conc_file'
    }

    # else:
    #   for i in range(len(query_list)):
    #     batch_proteins = query_list[i:i+500]
    #     hippie_query_api(batch_proteins, confidence_score)


    response = requests.post(url,params=parameters) #make a query in HIPPIE using the uniprot ids from the file
    print ('\n Response:', response.status_code)
    if response.status_code == 200:

      soup = BeautifulSoup(response.content, 'html.parser')# Parse the HTML response to find the form
      form = soup.find('form', {'id': 'netQuery'})

      if form:
        # Extract form action URL and hidden input fields
        action_url = form.get('action')
        form_data = {input_tag.get('name'): input_tag.get('value') for input_tag in form.find_all('input')}
        # print ('\n form_data:', form_data)
        # print ('\n action_url:', action_url)
        full_action_url = requests.compat.urljoin(response.url, action_url)

        print(f"Submitting the next form to: {full_action_url}")

    # Step 4: Submit the form to download the file
        download_response = requests.post(full_action_url, data=form_data, stream=True)
        print ('Download HIPPIE results: ', response.url) #this is a link that downloads the tsv file when it is clicked - i suppose that it is the same link that they sent to us via email
        print ('Download_response:', download_response)

      else: #if form not found
        download_response = "No form found on the page."
        print("\n No form found on the page.")

    elif str(response.status_code) in error_messages.keys() :
      print ('\n %s: %s' % (response.status_code, error_messages[str(response.status_code)]))

    else: #if status code != 200 and there is not such response in error messages
      print ('\n Error:', response.status_code)

    # filename = 'hippie_results.tsv'
    #if I create a csv file to save my results, again it saves the query uniprot ids in an HTML format
    # urlretrieve(new_url, filename)
    with open(filename, 'w') as file:
        file.write(download_response.content.decode('utf-8'))


    # df = pd.DataFrame(download_response.content.decode('utf-8'))
    # # df.to_csv(filename, index=False, sep='\t')
    # print (df)


    print (f'\n The file {filename} has been saved!')

    # read_hippie_results(filename)

def read_hippie_results (input_filename, output_filename):

  header = ['uniprot id 1', 'entrez gene id 1', 'gene name 1', 'uniprot id 2', 'entrez gene id 2', 'gene name 2', 'score'] #these are the columns of the conc file that hippie db returns
  index = 0
  #in case that a protein is not found the message "didn't find (.*) in the database, skipped it" appears at the top of the file before the results that start from header. So with the code below I try to obtain both the proteins not found and the results from HIPPIE

  with open(input_filename, 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    # print (list(reader))
    for i, row in enumerate(reader):
      if row != header:
        if row != []:
          protein = ((row[0]))# if i write print (row) it will print a list, so in this way it prints a str with a message with the protein which was not found
          pattern = r"didn't find (.*) in the database, skipped it"

          matches = re.finditer(pattern, protein)
          for m in matches:
              protein = (m.group(1)) #m.group is a string

          with open ('not_found_proteins.txt', 'a') as file:
            if protein not in 'not_found_proteins.txt':
              file.write(f'%s\n' % (protein))
      else:
        index = i
        # print (index)
        break

  hippie_interactions_df = pd.read_table(input_filename, sep ='\t', skiprows = index)
  print (hippie_interactions_df.head().to_string())
  print (hippie_interactions_df.columns.values)


  hippie_interactions_df.dropna(inplace = True, subset = ['uniprot id 1', 'uniprot id 2']) #remove any rows with 'nan' values - note: 'nan' is considered as float

  #Prepare the columns for the creation of a new csv file
  confidence_values_list = list(hippie_interactions_df.loc[:, 'score'])
  interactor_A_list_uniprot_ids = list(hippie_interactions_df.loc[:, 'uniprot id 1']) #discuss why i have chosen
  interactor_B_list_uniprot_ids = list(hippie_interactions_df.loc[:, 'uniprot id 2'])

  print ('\n Interactor A uniprot ids:', len (interactor_A_list_uniprot_ids))
  print (interactor_A_list_uniprot_ids)
  print ('\n Interactor B uniprot ids:', len (interactor_B_list_uniprot_ids))
  print (interactor_B_list_uniprot_ids)

  hippie_protein_pairs = []

  for a, b in zip (interactor_A_list_uniprot_ids, interactor_B_list_uniprot_ids):
    pair = (a, b) #just a note here: if i remove the pair (b,a) which is the same with (a,b) maybe i lose biderectional associations - useful (?) for networks
    if pair not in hippie_protein_pairs:
      hippie_protein_pairs.append(pair)

  #then i convert the uniprot ids to the equivalent gene names - just a note: The uniprot ids here are in the form of 'B2L11_HUMAN' and not accession numbers
  # interactor_A_list_genes = (protein_to_gene_conversion(interactor_A_list_uniprot_ids)['Gene'])

  dict_genes = protein_to_protein_accession_conversion(interactor_A_list_uniprot_ids + interactor_B_list_uniprot_ids)
  # interactor_B_dict_genes = uniprot_id_mapping(interactor_B_list_uniprot_ids)

  hippie_gene_pairs = []
  scores = []
  interactor_A_list_genes = []
  interactor_B_list_genes = []

  for pair in hippie_protein_pairs:
    protein_a = pair[0]
    protein_b = pair[1]
    index = hippie_protein_pairs.index(pair)

    if protein_a in dict_genes.keys() and protein_b in dict_genes.keys():
      gene_a = dict_genes[protein_a]['Acc']
      gene_b = dict_genes[protein_b]['Acc']

      if gene_a != 'None' and gene_b != 'None':
        pair = (gene_a, gene_b)

        if pair not in hippie_gene_pairs:
          hippie_gene_pairs.append(pair)
          interactor_A_list_genes.append(gene_a)
          interactor_B_list_genes.append(gene_b)
          scores.append(confidence_values_list[index])


  print ('\n Hippie gene pairs: ', len(hippie_gene_pairs))
  print ('\n Scores ', len(scores), '\n', scores)


  #create a new file that contains the interactors and the confidence value
  column_labels = ['Interactor A', 'Interactor B', 'Confidence Value']
  column_data = [interactor_A_list_genes, interactor_B_list_genes, scores]
  filename = output_filename
  new_csv_file(column_labels, column_data, filename, True)

  return hippie_gene_pairs #this is a list with tuples, each one represents a pair


def new_ppi_network(nodes, edges, filename):
    pairs_list = []
    for pair in edges:
        # print (pair)
      if (pair[0] in nodes and pair[1] in nodes) and (pair[0] != 'None' or pair[1] != 'None'):
        query_protein_index = nodes.index(pair[0])  # pair is a tuple
        reference_protein_index = nodes.index(pair[1])  # pair is a tuple; the

        pairs_list.append((query_protein_index, reference_protein_index))
    dt = np.dtype('int', 'int')
    new_matrix = np.array(pairs_list, dtype=dt)

    # print(pairs_list)

    g = ig.Graph(pairs_list)
    # print(g)

    g.vs['nodes'] = nodes
    # g.vs['localization'] = localization
    # g.es['protein'] = edges
    layout = g.layout("kamada_kawai", maxiter = 1000)
    # col_dict = {
    #     'cytosol': 'blue',
    #     'nucleus': 'green',
    #     'mitochondrion': 'yellow',
    #     'None': 'white'
    # }
    # g.vs['color'] = [col_dict[loc] for loc in g.vs['localization']]

    fig, ax = plt.subplots()
    ig.plot(g, layout=layout, vertex_label =g.vs['nodes'], vertex_size = 35, target=ax, bbox=(1000, 2000), margin=2250)

    # plt.figure(figsize=(10, 6))

    plt.savefig(filename)
    plt.show()

    print(f'\n The network {filename} has been created!')

with open ('mobidb_not_found.txt', 'w') as file:
    pass

def mobidb_query (query_list, input_type):

  mobidb_query.counter +=1

  url = 'https://mobidb.org/api/download_page?'

  error_messages = {
      '400':' Bad Request: The request was invalid.',
      '401': 'Unauthorized: Authentication failed.',
      '404': 'Not Found: The requested resource does not exist.',
      '414': 'Request-URI Too Long: The request URI is too long.',
      '500': 'Internal Server Error: The server encountered an error.'
  }

  query = ','.join(query_list)

  if input_type == 'genes':
    parameters = {
        'gene': query,
        'ncbi_taxon_id': '9606',
        'prediction-disorder-alphafold': 'exists'
    }

  elif input_type =='proteins':
    parameters = {
        'acc': query,
        'ncbi_taxon_id': '9606',
        'prediction-disorder-alphafold': 'exists',
    }

  response = requests.get (url, parameters)
  print ('\n Response:', response)
  results_dict = response.json() #this is a dictionary
  print ('\n Response json: ', response.json())


  entries_dictionary = {}
  uniprot_id_list = []

  if response.status_code == 200:

    for entry in (results_dict['data']): #results_dict['data'] is a list, entry is a dict

      # entry = response_dict['data'][i] #this is a dictionary
      gene_name = entry['gene']
      uniprot_id = entry['acc']
      protein_length = entry['length']

      if 'reviewed' in entry.keys():
      #   if entry['reviewed'] == True:
        uniprot_id_list.append(uniprot_id)

        if 'prediction-disorder-alphafold' in entry.keys():
          # disprot = response_dict_data['curated-disorder-disprot'] #from the dictionary of dictionaries i select only the dictionary that is realated to disprot. the disprot dictionary is a dictionary of dictionaries
          alphafold_idr_content = entry['prediction-disorder-alphafold']['content_fraction'] # this is the value of the content_fraction dictionary and contains the idr content
          # print (response_dict_data['localization'])

          #I noticed that some genes have also synonyms, so if i take the gene_name as it is in some cases i also take the synonyms (e.g. 'VCPSynonyms=HEL-220' )

          patterns = [r'(.*)Synonyms(.*)', r'(.*)Name(.*)']

          for pattern in patterns:

            match = re.search(pattern, gene_name)  # Using `search` instead of `findall`
            if match:
              gene_name = match.group(1)

          entries_dictionary[uniprot_id] = {'Gene': gene_name, 'Alphafold IDR Content':float(f'{(alphafold_idr_content*100):2f}'), 'Length': protein_length}

        else:
          entries_dictionary[uniprot_id] = {'Gene': gene_name, 'Alphafold IDR Content': 0.0, 'Length': protein_length} # it means that is not available/not found yet, so i did not put 0.0

        if 'localization' in entry.keys():
          localization = entry['localization']
          entries_dictionary[uniprot_id]['Localization'] = localization
        else:
          entries_dictionary[uniprot_id]['Localization'] = 'None'
        # with ('mobidb_not_found.txt', 'a') as file:
        #   file.write(f'{uniprot_id} \n')
        #   # print ('\n', uniprot_id)

      else:
        print (uniprot_id, '\n')
        with open ('mobidb_not_found.txt', 'a') as txtfile:
          txtfile.write(f'{uniprot_id} \n')
        # with open (f'no_reviewed_proteins_mobidb_({mobidb_query.counter}).txt', 'a') as txtfile:
        #   txtfile.write(f'{uniprot_id} \n')


  elif str(response.status_code) in error_messages.keys() :
    print ('\n %s: %s'% (response.status_code, error_messages[str(response.status_code)]))


  else:
    print ('\n Error:', response.status_code)

  print ('\n Mobidb Results:', entries_dictionary)

  # df = pd.DataFrame(entries_dictionary)
  # df.to_csv('mobidb_results.csv')

  return entries_dictionary

def protein_interactors_counter (interactors_list, pairs_list):

  '''
  This function takes as argument a list with protein pairs as tuples and returns a dictionary
  that has keys the protein ids and as values a list with proteins that each protein interacts with
  '''

  # counts_dict = {}  #this is a dict that has as keys the gene name and as values the number of interactors
  interactors_dict = {} #it has as key the protein id and as value a list with the protein ids of the proteins with which interacts
  # pair_genes = []
  # for pair in pairs_list:
  #   interactor_a = pair[0]
  #   interactor_b = pair[1]

  #   if interactor_a not in pair_genes:
  #     pair_genes.append(interactor_a)
  #   if interactor_b not in pair_genes:
  #     pair_genes.append(interactor_b)

  # converted_genes_dict = gene_to_protein_conversion(pair_genes)
  # print (converted_genes_dict)

  for pair in pairs_list:
    interactor_a = pair[0]
    interactor_b = pair[1]

    # interactor_a = converted_genes_dict[gene_a]['Protein']
    # interactor_b = converted_genes_dict[gene_b]['Protein']

    if interactor_a in interactors_list: # since the interactors list are the query proteins that i used as my initial datasets I know that there are no 'None' values so is unesessary to check for 'None' values
      if interactor_a not in interactors_dict.keys():
        interactors_dict[interactor_a] = [interactor_b]
      elif interactor_b not in interactors_dict[interactor_a]:
        interactors_dict[interactor_a].append(interactor_b)

    if interactor_b in interactors_list:
      if interactor_b not in interactors_dict.keys():
        interactors_dict[interactor_b] = [interactor_a]
      elif interactor_a not in interactors_dict[interactor_b]:
        interactors_dict[interactor_b].append(interactor_a)

    # if interactor_a != 'None' or interactor_b != 'None':
    #   if interactor_a in interactors_list and interactor_a not in interactors_dict.keys(): # i create a dict with keys only the proteins that i used in my queries,
    #   #and not any other proteins that were found in the intermediate
    #     interactors_dict[interactor_a] = [interactor_b]
    #   elif interactor_b not in interactors_dict[interactor_a]:
    #     interactors_dict[interactor_a].append(interactor_b)

    #   if interactor_b in interactors_list and interactor_b not in interactors_dict.keys():
    #     interactors_dict[interactor_b] = [interactor_a] #maybe it counts itself
    #   elif interactor_a not in interactors_dict[interactor_b]:
    #     interactors_dict[interactor_b].append(interactor_a)
    # else:


  no_interactors = list (set(interactors_list).difference(set(interactors_dict.keys()))) #these are proteins from the interactors list that were found not to have any pair
  for interactor in no_interactors:
    interactors_dict[interactor] = []


  print ('\n Total interactors: ', len(interactors_dict.keys()))
  print ('\n Proteins with no interactors: ', len(no_interactors))

  #Save my results in files
  protein_interactors_counter.counter += 1
  filename = f'not_interacting_proteins_{protein_interactors_counter.counter}.txt'
  with open (filename, 'w') as file:
    for protein in no_interactors:
      file.write(f'{protein} \n')

  print (f'\n The file {filename} has been saved!')

  new_csv_file (column_labels = ['Interactor A', 'Interactor B'],
                column_data = [interactors_dict.keys(), interactors_dict.values()],
                filename = f'interacting_proteins_{protein_interactors_counter.counter}.csv',
                index_status = True)

  return interactors_dict


def read_multiunired_results (filename, score, output_type, reverse):
  '''
    tsv_file: a table with the results from multiunired
    score: takes float values of 1.0, 0.5 or 0.0 indicating the presence of an association, paralogues, or no interaction
    reverse:
  '''
  read_multiunired_results.counter += 1 #track the number of times the function is called

  if re.search(r'\.csv$', filename): #check if the file ends with '.csv', so the $ symbol is used to define that this pattern is at the end of the given string
    df = pd.read_csv(filename, index_col = 'Unnamed: 0')

  elif re.search(r'\.tsv$', filename):
    df = pd.read_csv(filename, sep='\t', index_col = 'Unnamed: 0')

  elif re.search(r'\.xlsx$',filename):
    df = pd.read_excel(filename, index_col = 'Unnamed: 0')



  df = df.drop('Overall_Score', axis = 1, inplace = False) # axis 1 are the columns, inplace means to change the existing df, so since i need only the protein names i remove the last column that contains the total score
  if reverse == True:
    df = df.T
  #this file is a table where each row represents the query protein and each column the reference protein

  rows = list(df.index)  # saves the proteins in a list (the rows in MultiUnired analysis contain the proteins from DisProt)
  columns = list(df.columns.values)

  print ('\n MultiUniRed Rows: ', len(rows))
  print ('\n MultiUniRed Columns: ', len(columns))

  #MultiUniReD rows and columns contain the proteins in this form 'Q68D10 (SPTY2D1)', so i use regex to take only the uniprot accession
  multiunired_rows = []
  multiunired_columns = []
  pattern = r'(.*) (.*)'

  if output_type == 'genes':
    index = 2
  elif output_type == 'proteins':
    index = 1

  for row in rows:
    match = re.search(pattern, row)
    if match:
      row = match.group(index) #if row=match.group(1) i take the protein name elif row=match.group(2) i take the gene name
    multiunired_rows.append(row)


  for column in columns:
    match = re.search(pattern, column)
    if match:
      column = match.group(index)
    multiunired_columns.append(column)


  scores_array = df.values.astype(float) #take all the values of the table
  interactions_matrix = np.nonzero(scores_array == score) #returns a tuple arrays, one for each dimension,
  # with the indices of the values that to the given condition (score)

  list_of_coordinates = list(zip(interactions_matrix[0],interactions_matrix[1])) # create a list of the coordinates by merging the two arrays--> we have a list of tuples
  #interacted[0] are the query proteins (rows), interacted[1] are the reference proteins (columns)

  total_multiunired_pairs = [] # a list to save the protein pairs
  protein_interactors_list = [] # a list to save all the interactors - pre-process for protein-to-gene conversion
  interactor_a_list = []
  interactor_b_list = []

  for coord in list_of_coordinates: #coord is a tuple
      # print (coord)
      query_protein_index = int(coord[0]) # coord[0] is a class 'numpy.int64' so i transform it into an int
      reference_protein_index = int(coord[1])

      #i take the index of the matrix that the score is equal to one, and i
      #use this index to find the gene at this position in the rows_genes list
      query_protein = multiunired_rows[query_protein_index]
      reference_protein = multiunired_columns[reference_protein_index]
      interactor_a_list.append(query_protein)
      interactor_b_list.append(reference_protein)

      pair = (query_protein, reference_protein)

      if query_protein not in protein_interactors_list:
        protein_interactors_list.append(query_protein)

      if reference_protein not in protein_interactors_list:
        protein_interactors_list.append(reference_protein)

      if pair not in total_multiunired_pairs:
          total_multiunired_pairs.append(pair)

  print ('\n Total protein pairs: ', len(total_multiunired_pairs))
  print (total_multiunired_pairs)

  print ('\n Total protein interactors:', len(protein_interactors_list))
  print (protein_interactors_list)

  new_csv_file(column_labels = ['Interactor A', 'Interactor B'], column_data = [interactor_a_list, interactor_b_list], filename = f'multiunired_interactions_{score}_({read_multiunired_results.counter}).csv', index_status = False)

  return total_multiunired_pairs

def pairs_intersection (pair_1, pair_2):
  '''
    pair_1, and pair_2 are lists of tuples that denote protein pairs
  '''
  #select the smaller list to reverse
  if len(pair_1) <= len(pair_2):
    sorted_list = pair_1
    unsorted_list = pair_2
  
  else:
    sorted_list = pair_2
    unsorted_list = pair_1
  
  #create the reverse list
  reversed_list = []
  
  for pair in sorted_list:
    reversed_tuple = tuple(reversed(pair))
    reversed_list.append(reversed_tuple)

  common = set(unsorted_list).intersection(sorted_list)
  print ('Common non-reversed pairs :', len(common))
  print (list(common))
  common_reversed = set(unsorted_list).intersection(reversed_list)
  print ('Common reversed pairs:', len(common_reversed))
  print (list(common_reversed))
  intersection = common.union(common_reversed)

  print ('Total common pairs:', len(intersection))
  print (list(intersection))
  
  duplicates = []
  

  return intersection
