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
from Bio import SeqIO
import re

# import annotation data
aroAnn = pd.read_csv('aro_categories_index.tsv', sep='\t')  # When reading, tsv files must have their delimiter stated
# privateAnn = pd.read_csv('private_models.csv')  # Read in annotations that are not used by CARD in order to remove
# them from the translated annotation list
aroIndex = pd.read_csv('aro_index.tsv', sep='\t')  # Read index file in order to compare gene  to gene family via
# protein accession
aroDB = list(SeqIO.parse("nucleotide_fasta_protein_homolog_model.fasta", "fasta"))


#//region Annotation Translation
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

def Diff(li1, li2): # get difference of two lists
    diff = (list(set(li1) - set(li2)))
    diff.sort()
    return diff

#//region Unsearchable annotation checking and Cut entries that can't be entered into the database
dupedRows = newAnn[newAnn.duplicated(subset=['DNA Accession', 'class', 'mechanism', 'group'], keep=False)].copy()  #
# Rows that are just duplicate annotations
multiRows = newAnn[newAnn.duplicated(subset=['DNA Accession', 'group'], keep=False)].copy()  # rows which have the same
# DNA accession and group
multiRows = multiRows[~multiRows['DNA Accession'].isin(dupedRows['DNA Accession'])]  # Removes entries which are
# complete duplicates of one another from multiRows, leaving only the ones that differ in everything except DNA
# Accession and group.If there are any entries in this dataframe, then searching with DNA Accession and group together
# will still result in multiple hits
protDupe = newAnn[newAnn.duplicated(subset=['Protein Accession'], keep=False)]  # find duplicate Prot. Acc.

# with pd.option_context('display.max_columns', 4):
#     print(multiRows)

# Provide the user with output detailing which entries were cut and why
cutEntries = list(multiRows.index)
cutLines = [x+2 for x in cutEntries]  # adds 2 to index to get line # in annotation source file (+1 from counting
# from 1 instead of 0, +1 from header)
print("\033[1;31;31m The following entries (line # in the annotation file) were cut from original file because their "
      "DNA accessions and groups overlapped, making it impossible to add them to the database file:\033[0;31;39m")
print(cutLines)
# print("(Debug) Their index values are: ")
# print(cutEntries)
newAnn.drop(cutEntries, axis=0, inplace=True)  # removes rows with duplicate groups and DNA Accessions

nullProt = list(newAnn[newAnn['Protein Accession'].isnull()].index)
nullProtLines = [x+2 for x in nullProt]
print("\033[1;31;31m The following entries (line # in the annotation file) were cut from original file because they "
      "lack a protein accession, making it impossible to add them to the database file:\033[0;31;39m")
print(nullProtLines)
# print("(Debug) Their index values are: ")
# print(nullProt)
newAnn.dropna(how='any', axis=0, subset=['Protein Accession'], inplace=True)  # removes rows with no protein
# accession



protDupeLines = [x+2 for x in list(protDupe.index)]
print("\033[1;31;31m The following entries (line # in the annotation file) were cut from original file because they "
      "have duplicate protein accessions, making it impossible to add them to the database file:\033[0;31;39m")
print(Diff(protDupeLines, nullProtLines))
# print("(Debug) Their index values are: ")
# print(Diff(list(protDupe.index),nullProt))  # only show those that were cut because of duplication, not null
newAnn.drop(protDupe, axis=0, inplace=True)  # removes rows with duplicate Protein Accessions

exit()

# TODO: Some entries have duplicate Prot. Acc. and different gene families (eg. 921/923), while others have duplicate
#  prot. Accessions, but identical gene families (and different DNA accessions). Cut all duplicate protein accessions
#  as PA is the only way to identify genes to gene families
newAnn.reset_index(drop=True, inplace=True)  # Reset index after dropping entries to prevent empty rows


# newAnn.to_csv('newAnn_test.csv')
# print(newAnn.loc[cutEntries])
# print(newAnn)


#//endregion

#//region Add gene from index file
indexCols = [6, 5]  # Columns for protein accession and gene from index file, respectively
# categoriesCols = [0, 2, 1]  # Protein accession, gene family, and DNA Accession columns from annotation file,
# respectively

# Create dataframes from the categories and index files to line up each entry's gene with its gene family by
# protein accession
# compCategories = aroAnn[aroAnn.columns[categoriesCols]].copy()
compIndex = aroIndex[aroIndex.columns[indexCols]].copy()
#
# compCategories = compCategories.drop(cutEntries, axis=0)
# compCategories = compCategories.dropna(how='any', axis=0)  # Culls rows with no protein accession(or any empty cell)
# compCategories.reset_index(drop=True, inplace=True)  # cut the same entries from the gene comparison DataFrame as were
# # cut from the annotation to maintain consistency between the two. This allows them to be compared by index instead
# # of contents (protein accession)

# print(compCategories)
print(compIndex)


def dataframe_difference(df1, df2, doc=False, which=None):  # Compares 2 dataframes for their contents and outputs
    # the results. set doc to true to output a csv file containing the merged dataframe. pass which='both'to check
    # those that are in both one and the other.
    comparison_df = df1.merge(df2,
                              indicator=True,
                              how='outer',
                              on='Protein Accession'
                              )
    if which is None:
        diff_df = comparison_df[comparison_df['_merge'] != 'both']  # returns only the entries which differ
    else:
        diff_df = comparison_df[comparison_df['_merge'] == which]
    if doc:
        diff_df.to_csv('diff.csv')
    return diff_df


if dataframe_difference(newAnn, compIndex, doc=False, which='left_only').index.size > 0:  # Checks that there are
    # no entries that are only in the annotation file
    print("\033[1;31;31m ERROR: Annotation is larger than index. Ensure you have the most recent version of both \033["
          "0;31;39m")
    exit()
# in the output file, "right_only" means that it is only in the index. 'left_only' should be empty, as the index
# should have all of the entries that are in the annotation, but not vice versa. Having things in 'right_only' is
# therefore fine, but having any in left_only means that the index does not contain all entries, and something has
# gone wrong with your download.

print(len(newAnn))
geneMerge = dataframe_difference(newAnn, compIndex, doc=True, which='both')


with pd.option_context('display.max_columns', 10):
    print(geneMerge[geneMerge.duplicated(subset=['Protein Accession'], keep=False)])
# TODO: BUG: dataframe_difference has more rows than newAnn, which shouldn't be possible. Entries are getting
#  duplicated somehow

# TODO: Take the entries that align to both, concatenate their "gene" column onto the annotation entry with the
#  matching protein accession

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
# TODO: Compare compIndex to newAnn, match Prot. Acc. from Index and newAnn, and then concatenate matching gene
#  onto newAnn
#//endregion


# TODO: Move this to just before creating the annotation file, if possible (Check to see if any Database translator
#  code relies on this code)
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

print('Debug: exit early. Check TODO.')
exit()

#//endregion

#//region Database Translation

annotationAccessions = list(newAnn['DNA Accession'])
annotationGroup = list(newAnn['group'])

dbAccessions = []
dbGene = []
p = re.compile('(gb\|)')  # find the string 'gb|'
for i in range(0,len(aroDB)):  # Create list of DNA Accessions from database file
    dbAccessions.append(p.sub('', aroDB[i].id))  # remove 'gb|' from each entry's header
    dbAccessions[i] = dbAccessions[i][:dbAccessions[i].index('|')]  # Cuts everything after the first |, leaving only
    # the DNA accession
    # TODO: Get gene from DB's accession annotation file and assign using Protein accession
    dbGene.append(p.sub('', aroDB[i].id))  # remove 'gb|' from each entry's header
    splitGroup = dbGene[i].split("|")  # separates the header into individual sections
    dbGene[i] = splitGroup[4]  # adds the 'group' section of that header into the dbGene list

x = 1301  # test value that determines which annotation/databse entry pair to print
print(aroDB[x].id)  # print original id. If print(aroDB[x].id) matches its DNA Accession with print(newHeaders[x]),
# the translator is creating the list of translated headers in the right order

match = []
newHeaders = ['error'] * len(dbAccessions)  # Create list that will contain all translated headers in the correct order
# for the translated database. If one is not filled in, it will be listed as "error"

#//region -PRIME notation conversion
for i in range(0,len(dbGene)):  # DB file is in ' notation. Translated Annotation file is in -PRIME notation.
    # Converts dbGene to be in -PRIME annotation so that comparison will be accurate
    dbGene[i] = dbGene[i].replace('(?<! )\(', '')
    # replaces open paren only when within a word
    dbGene[i] = dbGene[i].replace('\'\'\)$', '-DPRIME')
    dbGene[i] = dbGene[i].replace('\'\)$', '-PRIME')  #
    # Replace ' or '' when next to a close paren and at the end of a line (single gene family)
    dbGene[i] = dbGene[i].replace('\'\'\);', '-DPRIME;')
    dbGene[i] = dbGene[i].replace('\'\);', '-PRIME;')
    # Replace ' or '' when next to a close paren and next to a semicolon (multiple gene families)
    dbGene[i] = dbGene[i].replace('(?<=\d)\)$', '')
    dbGene[i] = dbGene[i].replace('(?<=\d)\);', ';')
    # replaces ) when preceded by a number and (at the end of a string or semicolon-separated)
#//endregion

for i in range(0, len(annotationAccessions)):  # Sorts headers into same order as database by matching by DNA
    # accession and gene family
    for n in range(0, len(dbAccessions)):
        if annotationAccessions[i] == dbAccessions[n] and annotationGroup[i] == dbGene[n]:  # match by accession and
            # family
            match.append([i, n])  # give list of indices of accessions that match each other
            newHeaders[n] = newAnn['header'].loc[i]  # set list to contain all translated headers in the correct
            # order for the translated database

print(newHeaders[x])
# print(aroDB[x])
print("Matching Accessions: ")
print(match)

# Error checking before proceeding to the file writing stage
for i in range(0,len(newHeaders)):  # Checks that all database entries have been assigned a header
    if newHeaders[i] == 'error':
        print("Database entry " + str((i+1)*2 - 1) + " has not been given a value")  # Indicates the entries with no
        # header
        errorPresent = True

if errorPresent == True:
    exit()


#//endregion

#//region Write final files

# Write annotation file
today = date.today()  # get current date
filename = ("CARD_to_AMRplusplus_Annotation_" + today.strftime("%Y_%b_%d") + ".csv")
# Exports AMR++-ready annotation file and names it based on the present date

newAnn = newAnn['header', 'type', 'class', 'mechanism', 'group']
print(newAnn)  # TODO: Does this output an annotation file with the proper columns?

print("EXIT Early. check TODO")
exit()
pd.DataFrame.to_csv(newAnn, filename, index=False)  # exports converted annotation file as csv

# Write Database file
newAroDB = SeqIO.parse("nucleotide_fasta_protein_homolog_model.fasta", "fasta")
translatedFilename = ("./CARD_to_AMRplusplus_Database_" + today.strftime("%Y_%b_%d") + ".fasta")

i = 0
with open(translatedFilename, 'w') as translated:
    for record in newAroDB:
        print("old:" + record.id)
        record.id = newHeaders[i]  # Changes header to translated header
        # record.description = newHeaders[i]
        record.description = ''
        print("NEW:" + record.id)
        i += 1
        SeqIO.write(record, translated, 'fasta-2line')  # writes fasta file line-by-line. 'fasta-2line' instead of
        # 'fasta' forces it to avoid using line breaks
# //endregion
