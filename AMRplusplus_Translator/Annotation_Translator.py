# Translates CARD Annotation file into one that AMR++ can read

# TODO: Remove private annotations from translated. Compare private list to translated list by DNA Accession and remove
#  those that match
# TODO: Need to add ‘RequiresSNPConfirmation’ flags to the translator. Current ARO source files I’m using are homologue
#  model only. Translate variant model as well, add the flag, combine with translated annotation file.
# TODO: The fact that CARD sorts its drugs and gene families by semicolon separation in a single cell an issue will
#  cause false negatives. Fix.


import numpy as np
import pandas as pd
from datetime import date

# import annotation data
AroAnn = pd.read_csv('aro_categories_index.tsv', sep='\t')  # When reading, tsv files must have their delimiter stated

typeCol = ['Drugs'] * len(AroAnn.index)  # Creates a list of the string 'Drugs' with as many values as the
# annotation file has. MEGARes has a type column, but ARO (mostly) only deals with drugs
typeAnn = pd.DataFrame(typeCol, columns=['type'])  # Turns that list into a Dataframe

# Create new table containing AMR++-relevant data
AroCols = [1, 3, 4, 2]  # Important columns from ARO data.
newAnn = AroAnn[AroAnn.columns[AroCols]].copy()  # Creates a Dataframe containing the Columns from the ARO annotation file
# ARO's columns (branches) are in a different order than MEGARes, so this reorders them.
# Also Add DNA accession as junk data to fill in header temporarily - fills the same spot as Meg_#

newAnn.columns = ['DNA Accession', 'class', 'mechanism', 'group']  # sets names of columns of new annotation

#//region -PRIME notation conversion
# Removes parentheses in gene family and converts ' to -PRIME and '' to -DPRIME. eg. AAC(6') = AAC6-PRIME
newAnn['group'] = newAnn['group'].str.replace('(?<! )\(', '')
# replaces open paren only when within a word
newAnn['group'] = newAnn['group'].str.replace('\'\'\)$', '-DPRIME')
newAnn['group'] = newAnn['group'].str.replace('\'\)$', '-PRIME')  #
# Replace ' or '' when next to a close paren and at the end of a line (single gene family)
newAnn['group'] = newAnn['group'].str.replace('\'\'\);', '-DPRIME;')
newAnn['group'] = newAnn['group'].str.replace('\'\);', '-PRIME;')
# Replace ' or '' when next to a close paren and next to a semicolon (multiple gene families)
newAnn['group'] = newAnn['group'].str.replace('(?<=\d)\)$', '')
newAnn['group'] = newAnn['group'].str.replace('(?<=\d)\);', ';')
# replaces ) when preceded by a number and (at the end of a string or semicolon-separated)
#//endregion
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

pd.DataFrame.to_csv(newAnn,filename, index=False)  # exports converted annotation file as csv