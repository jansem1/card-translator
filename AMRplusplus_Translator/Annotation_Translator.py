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
aroIndex = pd.read_csv('aro_index.tsv', sep='\t')  # Read index file in order to compare gene  to gene family via
# protein accession

# TODO: Add protein accession to newAnn for comparison, then cut it off before writing the annotation file

# Create new table containing AMR++-relevant data
aroCols = [0, 1, 3, 4, 2]  # Important columns from ARO data.
newAnn = aroAnn[aroAnn.columns[aroCols]].copy()  # Creates a Dataframe containing the Columns from the ARO annotation
# file. ARO's columns (branches) are in a different order than MEGARes, so this reorders them.
# Also Add DNA accession to allow the database to be searched for matching entries - fills the same spot as Meg_###

# TODO: Compare index and categories prot. acc. and gene families, then create a new file for the database translator
#  with the gene and gene family so that the database translator can search by DNA Accession and gene,
#  then replace the original header

newAnn.columns = ['Protein Accession', 'DNA Accession', 'class', 'mechanism', 'group']  # sets names of columns of new
# annotation file

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

#//region Unsearchable annotation checking and Cut entries that can't be entered into the database
dupedRows = newAnn[newAnn.duplicated(subset=['DNA Accession', 'class', 'mechanism', 'group'], keep=False)].copy()  #
# Rows that are just duplicate annotations

multiRows = newAnn[newAnn.duplicated(subset=['DNA Accession', 'group'], keep=False)].copy()  # rows which have the same
# DNA accession and group

multiRows = multiRows[~multiRows['DNA Accession'].isin(dupedRows['DNA Accession'])]  # Removes entries which are
# complete duplicates of one another from multiRows, leaving only the ones that differ in everything except DNA
# Accession and group.If there are any entries in this dataframe, then searching with DNA Accession and group together
# will still result in multiple hits

# with pd.option_context('display.max_columns', 4):
#     print(multiRows)

# Provide the user with output detailing which entries were cut and why
cutEntries = list(multiRows.index)
cutLines = [x+2 for x in cutEntries]  # adds 2 to index to get line # in annotation source file (+1 from counting
# from 1 instead of 0, +1 from header)
print("\033[1;31;31m The following entries (line # in the annotation file) were cut from original file because their "
      "DNA accessions and groups overlapped, making it impossible to add them to the database file:\033[0;31;39m")
print(cutLines)
print("(Debug) Their index values are: ")
print(cutEntries)
# null_rows=newAnn[newAnn.isnull().any(axis=0)].copy()
# print(null_rows)
# print(newAnn[null_rows])
newAnn = newAnn.drop(cutEntries, axis=0)  # removes rows that cannot be entered into the database
newAnn = newAnn.dropna(how='any', axis=0, subset=['Protein Accession'])  # removes rows with no protein accession,
# and can thus not be entered into the database
newAnn.reset_index(drop=True, inplace=True)  # Reset index after dropping entries to prevent empty rows
# TODO: Provide a list of entries cut due to a lack of protein accession

# newAnn.to_csv('newAnn_test.csv')
# print(newAnn.loc[cutEntries])
# print(newAnn)

print('Debug: exit early. Check TODO.')
exit()

#//endregion

#//region add gene from index file
indexCols = [6, 5]  # Columns for protein accession and gene from index file, respectively
categoriesCols = [0, 2, 1]  # Protein accession, gene family, and DNA Accession columns from annotation file,
# respectively

# Create dataframes of protein accessions and gene families
compCategories = aroAnn[aroAnn.columns[categoriesCols]].copy()
compCategories = compCategories.drop(cutEntries, axis=0)
compCategories = compCategories.dropna(how='any', axis=0)  # Culls rows with no protein accession(or any empty cell)
compCategories.reset_index(drop=True, inplace=True)  # cut the same entries from the gene comparison DataFrame as were
# cut from the annotation to maintain consistency between the two. This allows them to be compared by index instead
# of contents (protein accession)

compIndex = aroIndex[aroIndex.columns[indexCols]].copy()
print(compCategories)
print(compIndex)

def dataframe_difference(df1, df2, which=None):
    comparison_df = df1.merge(df2,  # merges the two dataframes
                              indicator=True,
                              how='outer',
                              on='Protein Accession'
                              )
    if which is None:
        diff_df = comparison_df[comparison_df['_merge'] != 'both']  # returns only the entries which differ
    else:
        diff_df = comparison_df[comparison_df['_merge'] == which]
    comparison_df.to_csv('diff.csv')
    return diff_df


print(dataframe_difference(compCategories, compIndex, which=None))  # pass which='both' to check those that are in
# both one and the other. in the output file, "right_only" means that it is only in the index
#TODO: Take the entries that align to both, concatenate their "gene" column onto the annotation entry with the
# matching protein accession

print('Debug: exit early. Check TODO.')
exit()


# matchGene = []
# match = []
# for i in range(0, len(compCategories)):
#     for n in range(0, len(compIndex)):
#         if compCategories['Protein Accession'].loc[i] == compIndex['Protein Accession'].loc[n]:
#             # matchGene.append(compIndex['Protein Accession'].loc[n])
#             match.append([i, n])
# print(match)
# TODO: Compare compIndex to compCategories, match gene from Index onto Categories, and then concatenate matching gene
#  onto the annotation file. Either cut that column off in the database after it's been used to search and then
#  overwrite the annotation file, or create a new file just for comparison
#//endregion

typeCol = ['Drugs'] * len(newAnn.index)  # Creates a list of the string 'Drugs' with as many values as the
# annotation file has. MEGARes has a type column, but ARO (mostly) only deals with drugs
typeAnn = pd.DataFrame(typeCol, columns=['type'])  # Turns that list into a Dataframe so it can be concatenated to
# the new Dataframe which will contain the translated data

appendCol = newAnn['DNA Accession'].map(str) + "|" + typeAnn['type'].map(str) + "|" + newAnn['class'].map(str) + "|" + \
            newAnn['mechanism'].map(str) + "|" + newAnn['group'].map(str)  # concatenate columns to make header by
# mapping them to strings  and putting a | between each term
newAnn = pd.concat([appendCol, typeAnn, newAnn['class'], newAnn['mechanism'], newAnn['group']], axis=1) #
# concatenates all columns that must be in the final annotation
newAnn.columns = ['header', 'type', 'class', 'mechanism', 'group']  # rename columns to match with those of AMR++
print(newAnn)


# TODO: If the 'class' section contains a semicolon (multiple drugs), change the section of the string between the
#  second and third |, as well as the corresponding class column entry, into "multi-drug resistance"

#//region Write final file
today = date.today()  # get current date
filename = ("CARD_to_AMRplusplus_Annotation_" + today.strftime("%Y_%b_%d") + ".csv")
# Exports AMR++-ready annotation file and names it based on the present date

# TODO Don't forget to uncomment me when ready to output annotation file
pd.DataFrame.to_csv(newAnn, filename, index=False)  # exports converted annotation file as csv
#//endregion