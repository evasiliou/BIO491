import pandas as pd
import csv
from gprofiler import GProfiler

# convert tsv to csv file
def tsv_to_csv(tsv_file, csv_filename):
    # reading given tsv file
    csv_file = pd.read_csv(tsv_file, sep='\t')
    # converting tsv to csv
    csv_file.to_csv(csv_filename, index=False)


def read_csv(csv_file, index):
    uniprot_ids = []
    with open(csv_file, 'r') as csvFile:
        csvReader = csv.reader(csvFile)
        next(csvReader)

        for row in csvReader:
            uniprot_ids.append(row[index])

    return uniprot_ids

def csv_to_df(csv_file):
    df = pd.read_csv(csv_file)
    return df

def write_txt(txt_file, data_list):
    with open(txt_file, 'w') as txtFile:
        for data in data_list:
            txtFile.write(data)
            txtFile.write('\n')

def protein_to_gene_conversion(dictionary):
    '''
    The protein_to_gene_conversion function takes as argument the variable dictionary and returns a list with the gene names
    The variable dictionary saves the dictionary that has been created after reading a file (e.g. overall_scores)
    '''

    gp = GProfiler(
        user_agent='ExampleTool',  # optional user agent
        return_dataframe=True  # return pandas dataframe or plain python structures
    )
    print("\n Total proteins before conversion: ", len(dictionary))

    # convert the uniprot ids to the equivalent gene accession numbers. the results are saved in a df
    query_proteins = list(dictionary.keys())

    converted_df = gp.convert(organism='hsapiens',
                              query=query_proteins,
                              target_namespace='ENTREZGENE_ACC')

    print(converted_df.head().to_string())

    incoming_proteins = list(converted_df.loc[:, 'incoming'])
    query_genes = list(converted_df.loc[:, 'name'])

    # with the block below I read the list that contains the keys of the provided dictionary and after finding the appropriate index I append at the value 'Gene' (which is a list) the equivalent gene name
    for protein in incoming_proteins:
        index = incoming_proteins.index(protein)
        if query_genes[index] not in dictionary[protein]['Gene']:
            dictionary[protein]['Gene'].append(query_genes[index])

    # with the code below I pick only the first occurence of each gene name
    # query_genes = []
    # for gene in converted_genes:
    #     if gene not in query_genes:
    #         query_genes.append(gene)
    #
    # for protein in query_proteins:
    #     dictionary[protein]['Gene'].append()

    print("\nTotal proteins after conversion: ", len(query_genes))

    return query_genes


def gene_enrichment(dictionary, filename):
    query_genes = protein_to_gene_conversion(dictionary)

    gp = GProfiler(
        user_agent='ExampleTool',  # optional user agent
        return_dataframe=True  # return pandas dataframe or plain python structures
    )

    query = gp.profile(organism='hsapiens',
                       query=query_genes,
                       sources=["GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "WP", "HPA", "CORUM", "HP"],
                       # sources=["GO:MF", "GO:CC", "GO:BP"],
                       no_evidences=False)
    #
    print("\n Total proteins after enrichment: ", len(query))
    print("\n Enrichment Results: \n", query.head().to_string())
    # # print (dictionary)
    # #print(list(query.columns.values))
    #
    go_terms = list((query.loc[:, "native"]))  # pathways with their unique accession number
    print('\n Total GO terms:', len(go_terms))
    name = list((query.loc[:, 'name']))  # a description of the pathway
    print('\n Total names:', len(name))
    p_value = list(query.loc[:, "p_value"])
    print('\n Total p-values:', len(p_value))
    intersections = list(
        query.loc[:, "intersections"])  # a list with the group of proteins that participate in the equivalent pathway
    print('\n Total intersections:', len(intersections))
    #print(intersections)  # prints a list of lists
    #
    # print("\n Total proteins after enrichment: ", len(dictionary))
    # # print (dictionary)
    #
    new_dataframe = pd.DataFrame({
        # "Uniprot IDs": uniprot_ids,
        # "Gene": genes,
        # "Score": scores,
        "p-value": p_value,
        # "Autophagy": autophagy,
        "GO terms": go_terms,
        'Pathway(name)': name,
        'Intersections': intersections
    })
    new_dataframe.to_csv(filename, index=True)

    print("\n The file", f'{filename}', "has been saved!\n -------------------")


def read_gene_enrichment_results(filename):
    '''
    The find_comm function takes as argument the name of the file that contains the results from
    the gene enrichment analysis and returns a dictionary.
    '''

    csv_file = filename
    df = pd.read_csv(csv_file, converters={'Intersections': pd.eval})  # the converter here is applied for the data under the 'Intersections' column
    # since they were saved as str and not as lists
    print(df.head().to_string())
    print(df.columns.values)
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


#Remove all the autophagy-related proteins from the DisProt dataset by comparing these proteins to the ones downloaded from UniProt
#convert the tsv file with scores into a csv file
tsv_file = 'disprot.tsv'
csv_file = 'disprot.csv'
#csv_file = 'disprot_output_100_sample.csv'

tsv_to_csv(tsv_file, csv_file)
disprot_df= csv_to_df(csv_file)

#tsv_to_csv(tsv_file, csv_file)
tsv_file = 'UniProt_proteins.tsv'
csv_file = 'UniProt_proteins.csv'
tsv_to_csv(tsv_file, csv_file)
autophagy_df= csv_to_df(csv_file)
# print (autophagy_df.head().to_string())
# print (len(list(autophagy_df.loc[:, 'Entry'])))

disprot_acc= []
for i in disprot_df.index:
    if disprot_df.loc[i,'acc'] not in disprot_acc: #create a list that contains all the idr proteins as in disprot DB with no duplicates
        disprot_acc.append(disprot_df.loc[i,'acc'])
print ('\n Total disprot proteins: ', len(disprot_acc))

autophagy_acc= []
for i in autophagy_df.index:
    # create a list that contains autophagy related proteins as in UniProt without duplicates
    if autophagy_df.loc[i, 'Entry'] not in autophagy_acc:
        autophagy_acc.append(autophagy_df.loc[i, 'Entry'])
print ('\n Total autophagy proteins in uniprot: ', len(autophagy_acc))


#compare the two lists in order to remove the autophagy uniprot ids and insert the uncommon uniprot ids in a separate list
disprot_no_autophagy = []
disprot_autophagy = []
for uniprot_id in disprot_acc:
    if uniprot_id not in autophagy_acc and uniprot_id not in disprot_no_autophagy:
        disprot_no_autophagy.append(uniprot_id)
    if uniprot_id in autophagy_acc and uniprot_id not in disprot_autophagy:
        disprot_autophagy.append(uniprot_id)
print ('\n Total disprot proteins no autophagy: ', len(disprot_no_autophagy))
print ('\n Total disprot proteins involved autophagy: ', len(disprot_autophagy))


#prepare the files for multiunired
disprot_no_autophagy_df = pd.DataFrame({'acc': disprot_no_autophagy})
disprot_no_autophagy_df.to_csv('disprot_no_autophagy.csv',index=False)

print("\n The file", 'disprot_no_autophagy.csv', "has been saved!\n -------------------")

#------------------------------------------

#Read the files for gene enrichment


#The file with the disprot accession numbers (all the proteins related to autophagy have been removed)
csv_file = 'disprot_no_autophagy.csv' #the new file has been saved in csv format
df = pd.read_csv(csv_file, sep = '\t')  # convert the csv file into a pandas dataframe
#print (list(df.columns.values))
uniprot_ids = df.loc[:, "acc"]  # save the uniprot ids in a Series pandas data type

overall_scores = {}
for i in range (len(uniprot_ids)):
    overall_scores[df.loc[i,"acc"]] = {"Gene": [], "p-value": [], "GO terms": [], 'Pathway(name)': [], "Autophagy": [False] } #create a dict with keys the uniprot id in each row and values the equivalent score in each rowfor

#protein_to_gene_conversion(overall_scores)
gene_enrichment(overall_scores,'disprot_no_autophagy_gene_enrichment.csv' )

#The file with the uniprot accession numbers of the not interacting proteins after the multiunired analysis (as a control)

# csv_file = 'proteins_not_interacting_disprot2.csv' #the new file has been saved in csv format
# df = pd.read_csv(csv_file, sep = '\t')  # convert the csv file into a pandas dataframe
# #print (list(df.columns.values))
# uniprot_ids = df.loc[:, "acc"]  # save the uniprot ids in a Series pandas data type
#
# overall_scores = {}
# for i in range (len(uniprot_ids)):
#     overall_scores[df.loc[i,"acc"]] = {"Gene": [], "p-value": [], "GO terms": [], 'Pathway(name)': [], "Autophagy": [False] } #create a dict with keys the uniprot id in each row and values the equivalent score in each rowfor
#
# protein_to_gene_conversion(overall_scores)
# gene_enrichment(overall_scores,'not_interacting2_gene_enrichment.csv' )

#The file with uniprot proteins related to autophagy (reference file)
csv_file = 'uniprot.csv'
df = pd.read_csv(csv_file)  # convert the csv file into a pandas dataframe
#print (list(df.columns.values))
uniprot_ids = df.loc[:, "acc"]  # save the uniprot ids in a Series pandas data type

overall_scores = {}
for i in range (len(uniprot_ids)):
    overall_scores[df.loc[i,"acc"]] = {"Gene": [], "p-value": [], "GO terms": [], 'Pathway(name)': [], "Autophagy": [False] } #create a dict with keys the uniprot id in each row and values the equivalent score in each rowfor

#protein_to_gene_conversion(overall_scores)
gene_enrichment(overall_scores, 'uniprot_gene_enrichment.csv')

#find the common pathways after gene enrichment

#the read_gene_enrichemnt_results returns a dict like this: dictionary[go_term] = {'p-value': p_value, 'Pathway(name)': pathway,'Intersections': intersections}

uniprot_dict = read_gene_enrichment_results('uniprot_gene_enrichment.csv')#this is a dictionary
disprot_dict = read_gene_enrichment_results('disprot_no_autophagy_gene_enrichment.csv')#this is a dictionary

print ('\n Uniprot dictionary: \n', uniprot_dict)
print ('\n Disprot dictionary: \n', disprot_dict)

#filter the two dictionaries in order to take only the significant go terms meaning only the ones with p_value<=0.05
uniprot_significant = []
disprot_significant = []

print ('\n Uniprot significant go terms:')
for go_term in uniprot_dict.keys():
    p_value = uniprot_dict[go_term]["p-value"]
    if p_value<=0.05:
        uniprot_significant.append(go_term)
        #print (go_term, '|', p_value)
print (len(uniprot_significant))

print ('\n Disprot significant go terms:')
for go_term in disprot_dict.keys():
    p_value = disprot_dict[go_term]["p-value"]
    if p_value<=0.05:
        disprot_significant.append(go_term)
        #print (go_term, '|', p_value)
print (len(disprot_significant))
#the following code identifies the common GO terms and filters the dictionaries to retain only those GO terms as keys
uniprot_go_terms = set (uniprot_significant)
disprot_go_terms = set (disprot_significant)
common_go_terms = list(uniprot_go_terms.intersection(disprot_go_terms))
print ('\n Common GO terms: ', len (common_go_terms))
uniprot_dict_keys = list(uniprot_dict.keys())
disprot_dict_keys = list(disprot_dict.keys())

for key in uniprot_dict_keys:
    if key not in common_go_terms:
        uniprot_dict.pop(key)

print ('\Total uniprot go terms:', len(uniprot_dict))

for key in disprot_dict_keys:
    if key not in common_go_terms:
        disprot_dict.pop(key)

print('\n Total disprot go terms:', len(disprot_dict))

#the block below saves in the common_p_values list the p_values of each

#the block below searches each list of genes line by line to find common genes

#Just a note here: if the go terms are not necessary we can remove the common_dict and write a simpler code block

common_dict = {} #create a dictionary that contains only the common elements between the two datasets
for go_term in common_go_terms:
    common_dict[go_term] = []

print (common_dict)
for go_term in common_go_terms:
    uniprot_genes = uniprot_dict[go_term]['Intersections'] #this is a list
    disprot_genes = disprot_dict[go_term]['Intersections'] #this is a list
    common_genes = set(uniprot_genes).intersection(set(disprot_genes)) #this is a set

    if common_genes != {}: # check if the result of the intersection is not empty
        common_dict[go_term] = list(common_genes)
    else:
        common_dict[go_term] = ['None']

    #common_dict[go_term]['Common Pathway'] = uniprot_dict[go_term]['Pathway(name)']

print ('\n Common genes:' , len (common_dict), '\n', common_dict)

#create the rows for the new file
common_genes_list = []
#common_pathway= []
for go_term in common_go_terms:
    common_genes_list.append(common_dict[go_term])
    #common_cellular_process.append(common_dict[go_term]['Common Process'])

new_dataframe = pd.DataFrame({
    #"Uniprot IDs": uniprot_ids,
    #"Gene": genes,
    #"Score": scores,
    "GO terms": common_go_terms,
    #"p-value": common_p_value,
    #'Pathway(name)': common_cellular_process,
    'Common Intersections': common_genes_list
    #"Autophagy": common_autophagy,
})

#each time i give a different name to the file
new_dataframe.to_csv('common_gene_enrichment.csv', index=True)

print("\n The file 'common_gene_enrichment.csv' has been saved!")

#With the code below we create a new dictionary that has as keys the enriched genes and as values the go terms (maybe the processes are added later)

enriched_genes = {}

for go_term, common_genes in common_dict.items ():
    for gene in common_genes:
        if gene != None and gene not in enriched_genes.keys():
            enriched_genes[gene] = [go_term]
        else:
            enriched_genes[gene].append(go_term)
# print (common_dict.items()) #returns a tuple
# print (type(common_dict.get('GO:0010035')))
# print (common_dict.get('GO:0010035'))
print (common_dict.keys())
# for go_term in common_dict.keys():
#     gene_list = common_dict[go_term]
#     #print (type(gene_list), '\n')
#     if gene_list != []:
#         for gene in gene_list:
#             print (gene, '\n')
#             if gene not in enriched_genes.keys():
#                 enriched_genes[gene] = [go_term]
#             else:
#                 enriched_genes[gene].append(go_term)

#print (enriched_genes)
print (enriched_genes.keys())

#Below i filter the enriched_genes to remove the ones that do not have at least one equivalent go term
# genes_list = list (enriched_genes.keys())
# for gene in genes_list:
#     if enriched_genes[gene] == []:
#         enriched_genes.pop(gene)

print('\n Enriched genes:', len(list(enriched_genes.keys())), '\n', list(enriched_genes.keys()))

#The code block below reads the file with the results from multiunired and saves them in a dict

tsv_file = 'disprot_output2.tsv'
df=pd.read_table(tsv_file, sep='\t')

uniprot_ids = df.loc[:, "Unnamed: 0"]  # save the uniprot ids in a Series pandas data type
scores = df.loc[:, "Overall_Score"]  # save the scores in a Series pandas data type

# Gene enrichment using GProfiler
overall_scores = {}

for i in range (len(uniprot_ids)):
   overall_scores[df.loc[i,"Unnamed: 0"]] = {"Score": float(df.loc[i,"Overall_Score"]), "Gene": [], 'Enriched': []} #create a dict with keys the uniprot id in each row and values the equivalent score in each row

print (overall_scores)

#The file with the uniprot accession numbers of the not interacting proteins after the multiunired analysis (as a control)

# csv_file = 'proteins_not_interacting_disprot2.csv' #the new file has been saved in csv format
# df = pd.read_csv(csv_file, sep = '\t')  # convert the csv file into a pandas dataframe

# #print (list(df.columns.values))
# uniprot_ids = df.loc[:, "acc"]  # save the uniprot ids in a Series pandas data type
#
# overall_scores = {}
# for i in range (len(uniprot_ids)):
#     overall_scores[df.loc[i,"acc"]] = {"Gene": [], "p-value": [], "GO terms": [], 'Pathway(name)': [], "Enriched": [] } #create a dict with keys the uniprot id in each row and values the equivalent score in each rowfor
#
# protein_to_gene_conversion(overall_scores)
# gene_enrichment(overall_scores,'not_interacting2_gene_enrichment.csv' )

#the code below converts the uniprot ids into the equivalent genes

protein_to_gene_conversion(overall_scores) #print the list with the equivalent genes
print (overall_scores)
#print (uniprot_ids)


#Find the genes that are enriched

for uniprot_id in uniprot_ids:
    genes_list = list(overall_scores[uniprot_id]["Gene"]) #this is a list
    for gene in genes_list:
        if gene in list(enriched_genes.keys()):
            overall_scores[uniprot_id]["Enriched"].append(True)
            #print (overall_scores[uniprot_id]["Enriched"])
        else:
            overall_scores[uniprot_id]['Enriched'].append(False)

#create a new csv file
genes = []
enriched = []

for uniprot_id in uniprot_ids:
    genes.append (overall_scores[uniprot_id]['Gene'])
    enriched.append (overall_scores[uniprot_id]['Enriched'])
new_dataframe = pd.DataFrame({
    "Uniprot IDs": uniprot_ids,
    "Gene": genes,
    #"Score": scores,
    #"GO terms": common_go_terms,
    #"p-value": common_p_value,
    #'Pathway(name)': common_cellular_process,
    #'Common Intersections': common_genes_list
    "Enriched": enriched,
})

#each time i give a different name to the file
new_dataframe.to_csv('multiunired_results_with_enrichment.csv', index=True)

print("\n The file has been saved!")
