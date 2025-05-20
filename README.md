# Exploring the interactions of intrinsically disordered and autophagy-related proteins: a bioinformatics approach

This pipeline initiates with the selection of the two datasets which undergo duplication removal and gene enrichment analysis. Following this, they are used as queries in the MultiUniRed text mining tool. A new query is submitted in HIPPIE database using the proteins from these datasets in order to retrieve the experimentally verified interactors of these proteins. Afterwards, the results from text mining analysis and the experimental data are compared in order to find the protein associations either from text mining or experiments or the common interactions between the two. In the following step, plots showing the relation between the number of interactors and the disorder percentage or the length of the protein, as well as box plots are generated. Finally, proteins with interesting characteristics are selected and their experimentally or text-mining derived associations are visualized using PPINs. The methodology consists of twelve main steps (Figure 5), as described below.

Note that for all the analyses, except for gene enrichment, the UniProt acc (e.g.‘P12931’) is used and in cases that other identifiers are used (e.g. gene names or UniProt IDs), they are converted using the ID mapping tools. This is to avoid cases of gene synonyms, and we wanted to have a uniform identifier for all the analyses. Moreover, we consider interactions and associations synonymous terms. 

Step 1: Datasets Selection  

For this analysis the two datasets come from data stored in DisProt (Aspromonte et al. 2023) and UniProt (Consortium 2023). For the first dataset all the human proteins that contain IDRs were retrieved from DisProt using as keyword the NCBI Taxon ID ‘9606’ which refers to Homo sapiens. For the second dataset all the human reviewed proteins that are involved in autophagy or are related to autophagy were retrieved from UniProt using as keyword the GO term ID ‘go:0006914’. The data was downloaded on the 18th of February 2025.

Step 2: Duplication Removal and UniProt ID Extraction  

This step involves two processes that are necessary for gene enrichment that will be applied in the next step. Firstly, the UniProt Accession Number (hearafter UniProt acc) is extracted from both datasets and possible duplicate entries are removed. Finally, the UniProt acc are converted to the equivalent gene name using the UniProtMapper tool. 

Step 3: Gene Enrichment Analysis

In this step, gene enrichment analysis is applied on both datasets in order to find overrepresented biological processes using the gprofiler package. Conversion of UniProt acc to the equivalent gene name is requested before applying gene enrichment. Statistically significant associations considered those with a p-value ≤ 0.05. Following this, the enriched go terms, from both datasets are compared aiming to find common go terms/pathways in both datasets, which then we consider as markers of the biological processes that we study (autophagy). 

Step 4: Submitting a query in the MultiUniRed text mining tool

In this step, MultiUniRed literature mining tool is used in order to find possible interactions between IDPs and autophagy proteins. This step is repeated two times aiming to obtain the maximum possible associations. The first time we consider as query file the dataset with the disordered proteins and as reference file the autophagy proteins when the second time the datasets are uploaded in MultiUniRed in reverse order. For the purposes of this study, MultiUniRed was run at the backend and not using the online available tool since the dataset was very large. 

Step 5: Retrieval of protein-protein associations from MultiUniRed 

From the previous step we have two tables with the scores for the associations of each protein pair. We select only the pairs that have score equal to 1, and after checking for duplicates, we save the protein pairs into a tuple for further analysis. The protein pairs are also saved in a .csv file. 

Step 6: Retrieval of protein-protein interaction pairs from HIPPIE 

In this step, we take advantage of the HIPPIE’s API tool for querying the database. Since the dataset is large, we split it into smaller datasets with 200 entries and the process is repeated until all the proteins to be submitted. We also set a threshold to obtain the protein associations with a score at least equal to 0.63. The results from each query are returned into .tsv files where the first two columns represent the query protein with its UniProt id and its gene name and the next two columns the protein that interacts with again with its UniProt id and its gene name. There is also another column with the score of each interaction. In cases where a protein has not been found, it appears at the top of the .tsv file before the table with the associations. The protein pairs are collected, and then they are converted to the equivalent UniProt acc using the UniProt ID mapper. The final protein pairs are saved into a tuple for further analysis as well as in a .csv file.

Step 7: Disorder content - Query in MobiDB 

The MobiDB API tool is used to make queries using the UniProt acc of the proteins from the two datasets in order to collect data about the gene name, the length of the protein, and the pLDDT score from AlphaFold (cell localization). This information is collected only for the reviewed proteins. The results are saved in a .csv file, and the proteins that are excluded from the analysis are saved in a different .txt file. 

Step 8: Counting the Number of Interactors 

We created a function to count the number of interactors for each protein from the initial datasets (disordered and autophagy). This function takes as arguments, a list with the proteins which we examine their interactors, and a list with all the protein pairs, each one represented in a tuple format. A dictionary is created with keys the proteins of interest, and values lists of their interactors as found from the protein pairs. The length of the list is then calculated and represents the number of interactors for each protein. This function was called twice for MultiUniRed and HIPPIE protein pairs separately. 
Step 9: Comparison between the protein pairs 

The pairs from MultiUnired and HIPPIE analyses are further analyzed in order to find the associations that are only found either in MultiUnired or HIPPIE, as well as the common associations between the two datasets. For the unique interactions in the two datasets the difference set method was used, while for the common associations we created a function independently. This function takes as arguments two lists with pairs in tuples and checks which of the two lists is the smallest. After selecting the smaller one, it is reversed and then the protein pairs of the largest list are compared against the normal and the reversed dataset. In this way we managed to eliminate any duplicate associations. 

Step 10: Plots generation 

All the plots were generated using the matplotlib package. A Venn diagram is created showing the common interactions found experimentally or from the text-mining tool. Scatter plots are created showing the number of interactors in relation to the disordered content or the length of the proteins for the MultiUniRed and HIPPIE interactions separately. Additionally, we applied linear regression to investigate whether our data follow a trend. We also split the dataset into proteins that are disordered and ordered. We considered ordered proteins all the proteins that have AlphaFold prediction score equal to zero. Thus, histograms are created showing the frequency and the range of the length of the disordered and ordered proteins separately. Box plots were also generated illustrating the average number of interactors between disordered and ordered proteins for the MultiUniReD and HIPPIE interactions separately.   

Step 11: Network construction 

A function is used to select one or more proteins of interest in order to create a .csv file with all the interactors and the associations between them. This csv file is submitted in Cytoscape to create the network. The function takes as arguments a list with the protein(s) of interest and a list with tuples that represent the protein pairs that one wants to examine if the requested proteins are involved in. It also takes as arguments the input type (protein or gene) in order for the necessary transformations to be executed, and the name of the .csv file that will be created at the end. If there are any associations between the protein(s) of interest in the given dataset, the created .csv file is imported in Cytoscape and the network is generated. Each node represents a protein, and the edges represent the interaction between the protein and the interactor. For the current study some proteins with specific characteristics were selected for the visualisation of their interactions within the HIPPIE dataset, MultiUniReD dataset or both, since all these datasets were quite big and led to the formation of the so-called ‘hairballs’ when the networks were created. The nodes overlapped with one another and this made the understanding of the interactions and the extraction of any possible patterns or significant associations very difficult or even impossible.

