#the hippie_query function can be a backup in case that the hippie_api does not work
def hippie_query (query_list,confidence_score):

    '''
    This function takes as argument a file that contains only the enriched genes (previously created). From this file takes the uniprot ids and after making the queryin hippie it finds the interactions of these proteins that were found experimentally.
    It also saves the results from hippie in a tsv file 'hippie_results_0.63.tsv' where '0.63' is the confidence score that was used for the analysis.
    '''


    hippie_interactions_df = pd.read_table('HIPPIE_all_db.tsv' ) #saves in a df all the interactions in hippie
    # print (hippie_interactions_df)

    hippie_query_pairs = [] #this is a list with tuples that has the pairs of all the hippie interactions
    interactor_A = []
    interactor_B = []
    confidence_values_list = []

    for i in range (len(hippie_interactions_df)):
      confidence_value = float(hippie_interactions_df.loc[i, 'Confidence Value'])
      # print (type(confidence_value))
      interactor_a = hippie_interactions_df.loc[i, 'Gene Name Interactor A']
      interactor_b = hippie_interactions_df.loc[i, 'Gene Name Interactor B']

      if confidence_value >= float(confidence_score) and (interactor_a in query_list or interactor_b in query_list):
        pair = (interactor_a, interactor_b)
        if pair not in hippie_query_pairs:
          hippie_query_pairs.append(pair)

      interactor_A.append(interactor_a)
      interactor_B.append(interactor_b)
      confidence_values_list.append(confidence_value)

    print ('\n Total HIPPIE interactions: ', len(hippie_query_pairs))
    # print (len(hippie_query_pairs))


    column_labels = ['Interactor A', 'Interactor B', 'Confidence Value']
    column_data = [interactor_A, interactor_B, confidence_values_list]
    filename = 'hippie_results.csv'
    new_csv_file(column_labels, column_data, filename, True)


    return hippie_query_pairs #returns a list of tuples that contains the results from the query
