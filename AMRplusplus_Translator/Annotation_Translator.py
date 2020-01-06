# Translates CARD Annotation file into one that AMR++ can read

import numpy as np
import pandas as pd
from datetime import date

# import annotation data
AroAnn = pd.read_csv('aro_categories_index.tsv', sep='\t') # When reading, tsv files must have their delimiter stated
MegaAnn = pd.read_csv('megares_modified_annotations_v2.00.csv') # MEGARes Annotation

# print(AroAnn)
# print(MegaAnn)

#TODO: Add "type" column
#TODO:Is the fact that CARD sorts its drugs by semicolon separation in a single cell an issue? - WILL IT CAUSE FALSE NEGATIVES? eg. AMR++ doesn't list an item as carbapenem resistant if it is labelled "cephalosporn;penam"?
#TODO: Add type column and find a way to label ARO terms accordingly (Drug, multi-compound resistant, biocide, metal). Just label all as Drugs?
#TODO: import CARD database (nucleotide_fasta_...), extract the header, add CMG into header
#TODO: Convert AAC(2') and other such groups to -PRIME notation

# Create new table containing AMR++-relevant data
AroCols = [1,3,4,2] # Important columns from ARO data.
newAnn = AroAnn[AroAnn.columns[AroCols]] # Creates a Dataframe containing the Columns from the ARO annotation file
# ARO's branches are in a different order than MEGARes, hence 1,3,4,2.
# Also Add DNA accession (second column) as junk data to fill in header temporarily - fills the same spot as Meg_#

newAnn.columns = ['DNA Accession', 'class', 'mechanism', 'group'] # rename columns to match with those of AMR++

appendCol = newAnn['DNA Accession'].map(str) + "|" + newAnn['class'].map(str) + "|" + newAnn['mechanism'].map(str) + "|" + newAnn['group'].map(str) # concatenate columns to make header
appColumns = [1,2,3] # Sets the columns from newAnn to be copied over to the final CSV
finalAnn = pd.concat([appendCol, newAnn[newAnn.columns[appColumns]]], axis=1)
finalAnn.columns = ['header', 'class', 'mechanism', 'group']
print(finalAnn)

today = date.today()
filename = ("CARD_to_AMRplusplus_Annotation_" + today.strftime("%Y_%b_%d") + ".csv") # Exports AMR++-ready annotation file and names it based on the present date


pd.DataFrame.to_csv(finalAnn,filename, index=False) # exports converted annotation file as csv