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
aroAnn = pd.read_csv('aro_categories_index.tsv', sep='\t')  # When reading, tsv files must have their delimiter stated
# privateAnn = pd.read_csv('private_models.csv')  # Read in annotations that are not used by CARD in order to remove
# them from the translated annotation list

# Create new table containing AMR++-relevant data
aroCols = [1, 3, 4, 2]  # Important columns from ARO data.
newAnn = aroAnn[aroAnn.columns[aroCols]].copy()  # Creates a Dataframe containing the Columns from the ARO annotation file
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

#//region Unsearchable annotation checking
dupedRows = newAnn[newAnn.duplicated(keep=False)].copy()  # Rows that are just duplicate annotations
multiRows = newAnn[newAnn.duplicated(subset=['DNA Accession', 'group'], keep=False)].copy()  # rows which have the same
# DNA accession and group
# print(dupedRows)
multiRows = multiRows[~multiRows['DNA Accession'].isin(dupedRows['DNA Accession'])]  # Removes entries which are
# complete duplicates of one another, leaving only the ones that differ in everything except DNA Accession and group.
# If there are any entries in this dataframe, then searching with DNA Accession and group together will still result in
# multiple hits
# with pd.option_context('display.max_columns', 4):
    # print(multiRows)  # Print entries that are duplicates of one another (line number in annotation
    # file is multirows' index number + 2)

# print(newAnn[~newAnn['DNA Accession'].isin(multiRows['DNA Accession'])])
# newAnn = newAnn[~newAnn['DNA Accession'].isin(multiRows['DNA Accession'])]
cutEntries = list(multiRows.index)
cutLines = [x+2 for x in cutEntries]
print("\033[1;31;31m The following entries (line # in the annotation file) were cut from original file because their "
      "DNA accessions and groups overlapped, making it impossible to add them to the database file:\033[0;31;39m")
print(cutLines)
print("Their index values are: ")
print(cutEntries)
newAnn = newAnn.drop(cutEntries,axis=0)  # removes rows that cannot be entered into the database
newAnn.reset_index(drop=True, inplace=True)  # Reset index after dropping entries

# print(newAnn.loc[cutEntries])
# print(newAnn)
#//endregion

typeCol = ['Drugs'] * len(aroAnn.index)  # Creates a list of the string 'Drugs' with as many values as the
# annotation file has. MEGARes has a type column, but ARO (mostly) only deals with drugs
typeAnn = pd.DataFrame(typeCol, columns=['type'])  # Turns that list into a Dataframe so it can be concatenated to
# the new Dataframe which will contain the translated data

appendCol = newAnn['DNA Accession'].map(str) + "|" + typeAnn['type'].map(str) + "|" + newAnn['class'].map(str) + "|" + \
            newAnn['mechanism'].map(str) + "|" + newAnn['group'].map(str)  # concatenate columns to make header by
# mapping them to strings  and putting a | between each term
appColumns = [1, 2, 3]  # Sets the columns from newAnn to be copied over to the final CSV
newAnn = pd.concat([appendCol, typeAnn, newAnn[newAnn.columns[appColumns]]], axis=1)
newAnn.columns = ['header', 'type', 'class', 'mechanism', 'group']  # rename columns to match with those of AMR++

# TODO: If the 'class' section contains a semicolon (multiple drugs), change the section of the string between the
#  second and third |, as well as the corresponding class column entry, into "multi-drug resistance"

# print(newAnn)

today = date.today()  # get current date
filename = ("CARD_to_AMRplusplus_Annotation_" + today.strftime("%Y_%b_%d") + ".csv")
# Exports AMR++-ready annotation file and names it based on the present date

# TODO Don't forget to uncomment me when ready to output annotation file
pd.DataFrame.to_csv(newAnn,filename, index=False)  # exports converted annotation file as csv