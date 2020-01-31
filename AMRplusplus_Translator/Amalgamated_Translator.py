# Translates CARD Annotation file into one that AMR++ can read

# TODO: Remove private annotations from translated. Compare private list to translated list by DNA Accession and remove
#  those that match
# TODO: Need to add ‘RequiresSNPConfirmation’ flags to the translator. Current ARO source files I’m using are homologue
#  model only. Translate variant model as well, add the flag, combine with translated annotation file.
# TODO: The fact that CARD sorts its drugs and gene families by semicolon separation in a single cell an issue will
#  cause false negatives. Fix.
# TODO: Create README file for program


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
# 0 = Protein Accession; 1 = DNA Accession; 2 = Gene family; 3 = Drug; 4 = Class
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

def Diff(li1, li2): # entries that are in list 1 but not list 2
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
# Duplicate DNA accession and group cutting
cutEntries = list(multiRows.index)
cutLines = [x+2 for x in cutEntries]  # adds 2 to index to get line # in annotation source file (+1 from counting
# from 1 instead of 0, +1 from header)
print("\033[1;31;31m The following entries (line # in the annotation file) were cut from original file because their "
      "DNA accessions and groups overlapped, making it impossible to add them to the database file:\033[0;31;39m")
print(cutLines)
# print("(Debug) Their index values are: ")
# print(cutEntries)
newAnn.drop(cutEntries, axis=0, inplace=True)  # removes rows with duplicate groups and DNA Accessions

# Duplicate Protein Accession cutting
protDupeEntries = list(protDupe.index)
protDupeLines = [x+2 for x in protDupeEntries]
print("\033[1;31;31m The following entries (line # in the annotation file) were cut from original file because they "
      "lack or have duplicate protein accessions, making it impossible to add them to the database file:\033[0;31;39m")
print(protDupeLines)
# print("(Debug) Their index values are: ")
# print(protDupeEntries)  # only show those that were cut because of duplication, not null
newAnn.drop(protDupeEntries, axis=0, inplace=True)  # removes rows with duplicate Protein Accessions

# N/A cutting
nullEntries = newAnn[newAnn.isna().any(axis=1)].index.tolist()  # finds indices of any entry with any n/a cell
nullLines = [x+2 for x in nullEntries]
print("\033[1;31;31m The following entries (line # in the annotation file) were cut from original file because they "
      "lack a protein accession, DNA Accession, drug class, mechanism, or gene family, making it impossible to add "
      "them to the database file:\033[0;31;39m")
print(Diff(nullLines, protDupeLines))
# print("(Debug) Their index values are: ")
# print(Diff(nullEntries, protDupeEntries))
newAnn.dropna(how='any', axis=0, inplace=True)  # removes rows with no protein accession
newAnn.reset_index(drop=True, inplace=True)  # Reset index after dropping entries to prevent empty rows

dropTotal = len(Diff(nullEntries, protDupeEntries)) + len(protDupeEntries) + len(cutEntries)
percentDrop = round(100 * dropTotal/len(aroAnn), 2)
print("\n Total number of entries dropped: " + str(dropTotal) + ", which is " + str(percentDrop) + "% of total entries "
                                                                                                 "\n")
# newAnn.to_csv('newAnn_test.csv')
# print(newAnn.loc[cutEntries])
# print(newAnn)

#//endregion

#//region Add gene from index file


def dataframe_merge(df1, df2, doc=False, which=None, on='Protein Accession'):  # Compares 2 dataframes for their
    # contents and outputs the results. set doc to true to output a csv file containing the merged dataframe. pass
    # which='both'to check those that are in both one and the other.
    comparison_df = df1.merge(df2,
                              indicator=True,
                              how='outer',
                              on=on,
                              )
    duplicates = comparison_df[comparison_df.duplicated(keep='first')].index
    comparison_df.drop(duplicates, axis=0, inplace=True)
    comparison_df.reset_index(drop=True, inplace=True)  # Reset index after dropping entries to prevent empty rows
    if which is None:
        diff_df = comparison_df[comparison_df['_merge'] != 'both']  # returns only the entries which differ
    else:
        diff_df = comparison_df[comparison_df['_merge'] == which]
    if doc:
        diff_df.to_csv('diff.csv')
    return diff_df


indexCols = [6, 4]  # Columns for protein accession and gene, respectively, from index file
# Gene column is pulled from Model Name, not ARO name, because Model Name is used to build database headers
# Create dataframe from the index file to line up each entry's gene with its gene family by protein accession
compIndex = aroIndex[aroIndex.columns[indexCols]].copy()

if dataframe_merge(newAnn, compIndex, doc=False, which='left_only').index.size > 0:  # Checks that there are
    # no entries that are only in the annotation file
    print("\033[1;31;31m ERROR: Annotation is larger than index. Ensure you have the most recent version of both \033["
          "0;31;39m")
    exit()
# in the output file, "right_only" means that it is only in the index. 'left_only' should be empty, as the index
# should have all of the entries that are in the annotation, but not vice versa. Having things in 'right_only' is
# therefore fine, but having any in left_only means that the index does not contain all entries, and something has
# gone wrong with your download.


mergeCheck = dataframe_merge(newAnn, compIndex, doc=False, which='both')
if mergeCheck[mergeCheck.duplicated(subset=['Protein Accession'], keep=False)].index.size > 0:  # check for duplicates
    # after merge. Merge function sometimes produces duplicate entries
    print(mergeCheck[mergeCheck.duplicated(subset=['Protein Accession'], keep=False)].index.size)
    with pd.option_context('display.max_columns', 10):
        print(mergeCheck[mergeCheck.duplicated(subset=['Protein Accession'], keep=False)])  # print entries
        # duplicated by merging index and newAnn
    print("number of duplicates:" + str(len(mergeCheck[mergeCheck.duplicated(subset=['Protein Accession'], keep='first')])))
    print("newAnn length: " + str(len(newAnn)))
    print("mergeCheck length: " + str(len(mergeCheck)))
    print("ERROR: Duplicate entries detected. Exiting translator")
    exit()

newAnn = dataframe_merge(newAnn, compIndex, doc=True, which='both')  # add gene to newAnn by protein accession
newAnn.drop(['_merge'], axis=1, inplace=True)
# newAnn.reset_index(drop=True, inplace=True)  # Reset index after dropping entries to prevent empty rows

# with pd.option_context('display.max_columns', 10):
    # print(newAnn)

#//endregion
#//endregion

#//region Database Translation

annotationAccessions = list(newAnn['DNA Accession'])  # this needs to be here because DNA Accession gets cut from
# newAnn in the next section
annotationGene = list(newAnn['Model Name'])
# TODO: FIND A WAY TO SEARCH FOR DNA Accessions in newAnn. newAnn.str.find() gives a series.
#  - Maybe an issue with DNA Accessions being culled? If so, would probably see more database entries not matching
# newAnn['DNA Accession'].str.find()
print('Early Break. FIND A WAY TO SEARCH FOR DNA ACCESSION TO FIGURE OUT WHY 13 ENTRIES AREN\'T BEING matched' 
      'properly')
exit()
#//region Create AMR++-compliant header
typeCol = ['Drugs'] * len(newAnn.index)  # Creates a list of the string 'Drugs' with as many values as the
# annotation file has. MEGARes has a type column, but ARO (mostly) only deals with drugs
typeAnn = pd.Series(typeCol)  # Turns that list into a Dataframe so it can be concatenated to
# the new Dataframe which will contain the translated data

headerCol = newAnn['DNA Accession'].map(str) + "|" + typeAnn.map(str) + "|" + newAnn['class'].map(str) + "|" + \
            newAnn['mechanism'].map(str) + "|" + newAnn['group'].map(str)  # concatenate columns to make header

newAnn = pd.concat([headerCol, typeAnn, newAnn['class'], newAnn['mechanism'], newAnn['group'],
                    newAnn['Protein Accession'], newAnn['Model Name']], axis=1)  # concatenates all columns that must be
# in the final annotation

newAnn.columns = ['header', 'type', 'class', 'mechanism', 'group', 'Protein Accession', 'Model Name']  # rename columns
# to match with those of AMR++. Protein Accession and Model Name will be cut before export
# TODO: If the 'class' section contains a semicolon (multiple drugs), change the section of the string between the
#  second and third |, as well as the corresponding class column entry, into "multi-drug resistance"
#//endregion

#//region dbGenes and dbAccession Processor
dbAccessions = []
dbGenes = []
p = re.compile('(gb\|)')  # find the string 'gb|'
for i in range(0,len(aroDB)):  # Create lists of DNA Accessions and genes from database headers
    dbAccessions.append(p.sub('', aroDB[i].id))  # remove 'gb|' from each entry's header
    dbAccessions[i] = dbAccessions[i][:dbAccessions[i].index('|')]  # Cuts everything after the first |, leaving only
    # the DNA accession
    # TODO: When getting database gene, use aroDB.description (not aroDB.id) and then cut everything after the first [
    dbGenes.append(aroDB[i].description)
    # TODO: MATCHING STOPS WORKING IF THIS IS CHANGED TO .description. Does work if .description ir
    splitGroup = dbGenes[i].split("|")  # separates the header into individual sections
    dbGenes[i] = splitGroup[5]  # adds the 'group' section of the database header into the dbGenes list
    dbGenes[i] = dbGenes[i][:dbGenes[i].index(' [')]  # cuts species information, leaving only the gene
#//endregion

#//region -PRIME notation conversion for dbGenes
for i in range(0,len(dbGenes)):  # DB file is in ' notation. Translated Annotation file is in -PRIME notation.
    # Converts dbGenes to be in -PRIME annotation so that comparison will be accurate
    dbGenes[i] = dbGenes[i].replace('(?<! )\(', '')
    # replaces open paren only when within a word
    dbGenes[i] = dbGenes[i].replace('\'\'\)$', '-DPRIME')
    dbGenes[i] = dbGenes[i].replace('\'\)$', '-PRIME')  #
    # Replace ' or '' when next to a close paren and at the end of a line (single gene family)
    dbGenes[i] = dbGenes[i].replace('\'\'\);', '-DPRIME;')
    dbGenes[i] = dbGenes[i].replace('\'\);', '-PRIME;')
    # Replace ' or '' when next to a close paren and next to a semicolon (multiple gene families)
    dbGenes[i] = dbGenes[i].replace('(?<=\d)\)$', '')
    dbGenes[i] = dbGenes[i].replace('(?<=\d)\);', ';')
    # replaces ) when preceded by a number and (at the end of a string or semicolon-separated)
#//endregion


match = []
newHeaders = ['error'] * len(dbAccessions)  # Create list that will contain all translated headers in the correct order
# for the translated database. If one is not filled in, it will be listed as "error"

for i in range(0, len(annotationAccessions)):  # Sorts headers into same order as database by matching by DNA
    # accession and gene family
    for n in range(0, len(dbAccessions)):
        if annotationAccessions[i] == dbAccessions[n] and annotationGene[i] == dbGenes[n]:  # match by accession and
            # gene
            match.append([i, n])  # give list of indices of accessions that match each other
            newHeaders[n] = newAnn['header'].loc[i]  # set list to contain all translated headers in the correct
            # order for the translated database

# x = 305  # test value that determines which annotation/databse entry pair to print
x = 451

print(aroDB[x].description)  # print original id. If print(aroDB[x].id) matches its DNA Accession with print(newHeaders[x]),
# the translator is creating the list of translated headers in the right order
print(newHeaders[x])
# print(aroDB[x])
print("Matching Accessions: ")
print(match)
# exit()
# TODO: Cull Database entries with identical DNA Accessions and genes
# Error checking before proceeding to the file writing stage
noValue = 0
for i in range(0, len(newHeaders)):  # Checks that all database entries have been assigned a header
    if newHeaders[i] == 'error':
        print("Database entry " + str((i+1)*2-1) + " has not been given a value")  # Indicates the entries with no
        # header
        print("DEBUG: index = " + str(i))
        print(aroDB[i].description)
        noValue += 1
        errorPresent = True
print("Number of missing entries: " + str(noValue))
if errorPresent == True:
    exit()
# TODO: ARO NAMES IS NOT THE SAME AS THE DATABASE HEADER'S GENE SOMETIMES. SWITCH OVER TO USING "MODEL NAME"? Talk to
#  Brian/Andrew first

#//endregion

newAnn = newAnn['header', 'class', 'mechanism', 'group']  # drop all columns that are unneeded for annotation
# file
print(newAnn)
# TODO: Does this output an annotation file with the proper columns? Drop PA and Model Name

print("EXIT Early. check TODO")
exit()





#//region Write final files

# Write annotation file
today = date.today()  # get current date
filename = ("CARD_to_AMRplusplus_Annotation_" + today.strftime("%Y_%b_%d") + ".csv")
# Exports AMR++-ready annotation file and names it based on the present date


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
