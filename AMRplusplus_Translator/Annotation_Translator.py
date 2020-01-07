# Translates CARD Annotation file into one that AMR++ can read

# TODO: Need to label for multi-compound-resistant determinants for multidrug resistance,
#  or only for multiple types (eg. Drugs and metals)?
# TODO: Is the fact that CARD sorts its drugs by semicolon separation in a single cell an issue?
#  WILL IT CAUSE FALSE NEGATIVES? eg. AMR++ doesn't list an item as carbapenem resistant if it is
#  labelled "cephalosporin;penam"?
# TODO: Add type column and find a way to label ARO terms accordingly (Drug, multi-compound resistant, biocide, metal).
#  Just label all as Drugs?
# TODO: import CARD database (nucleotide_fasta_...), extract the header, replace with new headers. Will need to find
#  each header via DNA accession and replace it
# TODO: Convert AAC(2') and other such groups to -PRIME notation
#  - Find and delete all parentheses in each cell of the "group" column (BEFORE it is concatenated into the header)
#   - for open parens, make sure it excludes open parens with a space in front of them (Eg. see line 808 of Jan 6 export)
#   - For close parens, make sure it's ONLY at the end of the string
#  - Find all '' in the group column. Replace with '-DPRIME' string. Use {Dataframe}.str.replace
#  - Find all ' in the group column. Replace with '-PRIME' string. Use {Dataframe}.str.replac


import numpy as np
import pandas as pd
from datetime import date

# import annotation data
AroAnn = pd.read_csv('aro_categories_index.tsv', sep='\t')  # When reading, tsv files must have their delimiter stated
# MegaAnn = pd.read_csv('megares_modified_annotations_v2.00.csv') # MEGARes Annotation

typeCol = ['Drugs'] * len(AroAnn.index)  # Creates a list of the string 'Drugs' with as many values as the
# annotation file has. MEGARes has a type column, but ARO (mostly) only deals with drugs
typeAnn = pd.DataFrame(typeCol, columns=['type'])  # Turns that list into a Dataframe
# print(typeAnn)

# Create new table containing AMR++-relevant data
AroCols = [1, 3, 4, 2]  # Important columns from ARO data.
newAnn = AroAnn[AroAnn.columns[AroCols]].copy()  # Creates a Dataframe containing the Columns from the ARO annotation file
# ARO's columns (branches) are in a different order than MEGARes, so this reorders them.
# Also Add DNA accession as junk data to fill in header temporarily - fills the same spot as Meg_#

newAnn.columns = ['DNA Accession', 'class', 'mechanism', 'group']  # sets names of columns of new annotation

# Removes parentheses in gene family and converts ' to -PRIME. eg. AAC(6') = AAC6-PRIME
newAnn['group'] = newAnn['group'].str.replace('(?<! )\(', '')  # !!DANGER!! This will replace any open paren, anywhere. See TO-DO
# '(?<! ) is regex for "unless there is a space in front of the value"
newAnn['group'] = newAnn['group'].str.replace('\'\'\)$', '-DPRIME')
newAnn['group'] = newAnn['group'].str.replace('\'\)$', '-PRIME')  #
# Replace ' or '' when next to a close paren and at the end of a line (single gene family)

newAnn['group'] = newAnn['group'].str.replace('\'\'\);', '-DPRIME')  # will only replace ' when it is next to a close paren,
newAnn['group'] = newAnn['group'].str.replace('\'\);', '-PRIME')
# Replace ' or '' when next to a close paren and next to a semicolon (multiple gene families)

newAnn['group'] = newAnn['group'].str.replace('\)$', '')  # replaces ) at the end of a string
# ($ is a regex that causes this to only apply at the end)
newAnn['group'] = newAnn['group'].str.replace('\);', '') # replaces ) when in multiple gene families
# (see line 12 of converted file)
#TODO is this all the possible cases?

appendCol = newAnn['DNA Accession'].map(str) + "|" + typeAnn['type'].map(str) + "|" + newAnn['class'].map(str) + "|" + \
            newAnn['mechanism'].map(str) + "|" + newAnn['group'].map(str)  # concatenate columns to make header by
# mapping them to strings  and putting a | between each term
appColumns = [1, 2, 3]  # Sets the columns from newAnn to be copied over to the final CSV
newAnn = pd.concat([appendCol, typeAnn, newAnn[newAnn.columns[appColumns]]], axis=1)
newAnn.columns = ['header', 'type', 'class', 'mechanism', 'group']  # rename columns to match with those of AMR++

print(newAnn)

today = date.today()  # get current date
filename = ("CARD_to_AMRplusplus_Annotation_" + today.strftime("%Y_%b_%d") + ".csv")
# Exports AMR++-ready annotation file and names it based on the present date

# pd.DataFrame.to_csv(newAnn,filename, index=False)  # exports converted annotation file as csv
