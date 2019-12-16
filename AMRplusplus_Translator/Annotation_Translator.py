# Translates CARD Annotation file into one that AMR++ can read

import pandas as pd
from datetime import date

# import annotation data
AroAnn = pd.read_csv('aro_categories_index.tsv', sep='\t', index_col="DNA Accession") # When reading, tsv files must have their delimiter stated
MegaAnn = pd.read_csv('megares_annotations_v2.00.csv') # MEGARes Annotation

# print(AroAnn)
# print(MegaAnn)
#TODO: import CARD database (nucleotide_fasta_...), extract the header, add CMG into header
#TODO: Add this header to final conversion file as index column

# Create new table containing AMR++-relevant data
AroCols = [2,3,1] # Important columns from ARO data.
# ARO's branches are in a different order than MEGARes, hence 2,3,1.
# Add DNA accession (second column) as junk data to fill in header temporarily
NewAro = AroAnn[AroAnn.columns[AroCols]] # Creates a Dataframe containing the Columns from the ARO annotation file

# rename columns to match with those of AMR++
NewAro.index.names = ['header']
NewAro.columns = ['class', 'mechanism', 'group']

# Exports AMR++-ready annotation file
print(NewAro)
today = date.today()
filename = ("CARD_to_AMRplusplus_Annotation_" + today.strftime("%Y_%b_%d") + ".csv")
print(filename)
# pd.DataFrame.to_csv(NewAro,filename) # exports converted annotation file as csv