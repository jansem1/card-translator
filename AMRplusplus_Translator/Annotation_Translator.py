# Translates CARD Annotation file into one that AMR++ can read

import pandas as pd
from datetime import date

# import annotation data
AroAnn = pd.read_csv('aro_categories_index.tsv', sep='\t', index_col="DNA Accession") # When reading, tsv files must have their delimiter stated
MegaAnn = pd.read_csv('megares_modified_annotations_v2.00.csv') # MEGARes Annotation

# print(AroAnn)
# print(MegaAnn)
#TODO:Is the fact that CARD sorts its drugs by semicolon separation in a single cell an issue?
#TODO: Add type column and find a way to label ARO terms accordingly (Drug, multi-compound resistant, biocide, metal). Just label all as Drugs?
#TODO: import CARD database (nucleotide_fasta_...), extract the header, add CMG into header
#TODO: Add this header to final conversion file as index column
#TODO: Convert AAC(2') and other such groups to -PRIME notation

# Create new table containing AMR++-relevant data
AroCols = [2,3,1] # Important columns from ARO data.
# ARO's branches are in a different order than MEGARes, hence 2,3,1.
# Add DNA accession (second column) as junk data to fill in header temporarily
NewAnn = AroAnn[AroAnn.columns[AroCols]] # Creates a Dataframe containing the Columns from the ARO annotation file

# rename columns to match with those of AMR++
NewAnn.index.names = ['header']
NewAnn.columns = ['class', 'mechanism', 'group']

for i in range (0,len(NewAnn.index)):
    print(NewAnn.iloc[i,0])
#     print(NewAnn.index[i])

# Exports AMR++-ready annotation file and names it based on the present date
# print(NewAnn) # Displays table to verify that the translation is occurring properly
today = date.today()
filename = ("CARD_to_AMRplusplus_Annotation_" + today.strftime("%Y_%b_%d") + ".csv")
# print(filename)

# pd.DataFrame.to_csv(NewAnn,filename) # exports converted annotation file as csv