import pandas as pd
from datetime import date

# import annotation data
AroAnn = pd.read_csv('aro_categories_index.tsv', sep='\t', index_col="DNA Accession") # When reading, tsv files must have their delimiter stated
MegaAnn = pd.read_csv('megares_annotations_v1.01.csv')

# print(AroAnn)
# print(MegaAnn)

# Create new table containing AMR++-relevant data
AroCols = [2,3,1] # Important columns from ARO data.
# ARO's branches are in a different order than MEGARes, hence 2,3,1.
# Add DNA accession as junk data to fill in header temporarily
NewAro = AroAnn[AroAnn.columns[AroCols]]

# rename columns to match with those of AMR++
NewAro.index.names = ['header']
NewAro.columns = ['class', 'mechanism', 'group']

# Exports AMR++-ready annotation file
print(NewAro)
today = date.today()
filename = ("CARD_to_AMRplusplus_" + today.strftime("%Y_%b_%d") + ".csv")
print(filename)
# pd.DataFrame.to_csv(NewAro,filename)