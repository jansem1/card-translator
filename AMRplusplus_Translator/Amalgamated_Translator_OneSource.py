# Translates CARD Annotation file into one that AMR++ can read

# TODO: Remove private annotations from translated. Compare private list to translated list by DNA Accession and remove
#  those that match
# TODO: remove "_combined" from default filenames once this version takes over from the original translator
# TODO: Make the script operable from command line, and make it possible to change the names of the input files
# TODO: The fact that CARD sorts its drugs and gene families by semicolon separation in a single cell an issue will
#  cause false negatives. Fix.
# TODO: Remove "multi-drug resistance" if multi-group entries cannot be removed or changed to "multi-group" somehow
# TODO: Create README file for program

import numpy as np
import pandas as pd
from datetime import date
from Bio import SeqIO
import re
import argparse

# Define default filenames
aroIndexFile = 'aro_index_combined.tsv'
aroDBFile = "nucleotide_fasta_protein_homolog_model_combined.fasta"

#//region Parse commandline arguments to allow the user to input different filenames
parser = argparse.ArgumentParser()
# parser.add_argument("-a", help='Allows user to define the filepath of the ARO annotation file')
parser.add_argument("-i", help='Allows user to define the filepath of the ARO index file')
parser.add_argument("-d", help='Allows user to define the filepath of the ARO database file')

args = parser.parse_args()
# if isinstance(args.a, str):
#     aroAnnFile = args.a
# print("Annotation file: " + aroAnnFile)
if isinstance(args.i, str):
    aroIndexFile = args.i
print("Index file: " + aroIndexFile)
if isinstance(args.d, str):
    aroDBFile = args.d
print("Database file: " + aroDBFile)
#//endregion

#//region define functions

def Diff(li1, li2):  # entries that are in list 1 but not list 2
    diff = (list(set(li1) - set(li2)))
    diff.sort()
    return diff


def dataframe_merge(df1, df2, doc=False, which=None, on='Protein Accession', ind=True, byIndex=False):  # Compares 2
    # dataframes for their contents and outputs the results. set doc to true to output a csv file containing the
    # merged dataframe. pass which='both'to check those that are in both one and the other.
    comparison_df = df1.merge(df2,
                              indicator=True,
                              how='outer',
                              on=on,
                              left_index=byIndex,
                              right_index=byIndex
                              )
    duplicates = comparison_df[comparison_df.duplicated(keep='first')].index  # finds duplicate entries created by
    # the merge
    comparison_df.drop(duplicates, axis=0, inplace=True)  # drops duplicate entries created by the merge
    comparison_df.reset_index(drop=True, inplace=True)  # Reset index after dropping entries to prevent empty rows
    if ind:
        if which is None:
            merged_df = comparison_df[comparison_df['_merge'] != 'both']  # returns only the entries which differ
        else:
            merged_df = comparison_df[comparison_df['_merge'] == which]
    else:
        merged_df = comparison_df[comparison_df['_merge'] == 'both']
        # merged_df.drop(['_merge'], axis=1, inplace=True)
    if doc:
        merged_df.to_csv('diff.csv')
    return merged_df


#//endregion

# Import files
# aroAnn = pd.read_csv(aroAnnFile, sep='\t') # import annotation data
# privateAnn = pd.read_csv('private_models.csv')  # Read in annotations that are not used by CARD in order to remove
# them from the translated annotation list
aroIndex = pd.read_csv(aroIndexFile, sep='\t')  # Read index file in order to compare gene  to gene family via
# protein accession
aroDB = list(SeqIO.parse(aroDBFile, "fasta"))  # Read in CARD Database file as a list of seqRecord objects
newAroDB = SeqIO.parse(aroDBFile, "fasta")  # Read in CARD database as a Seqrecord generator, instead of a list,
# so that it can be overwritten later

#//region Annotation Translation

# Create new DataFrame containing AMR++-relevant data
aroCols = ['Protein Accession', 'DNA Accession', 'Drug Class', 'Resistance Mechanism', 'AMR Gene Family',
           'Model Name']  # Important
# columns from ARO data.
newAnn = aroIndex[aroCols].copy()  # Creates a Dataframe containing the Columns from the ARO annotation
# file. ARO's columns (branches) are in a different order than MEGARes, so this reorders them.
# Also Add DNA accession to allow the database to be searched for matching entries - fills the same spot as Meg_###

newAnn.columns = ['Protein Accession', 'DNA Accession', 'class', 'mechanism', 'group', 'Model Name']  # sets names of columns of new
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

#//region Unsearchable annotation checking and culling
dupedRows = newAnn[newAnn.duplicated(subset=['DNA Accession', 'class', 'mechanism', 'group'], keep=False)].copy()  #
# Rows that are just duplicate annotations
overlapRows = newAnn[newAnn.duplicated(subset=['DNA Accession', 'group'], keep=False)].copy()  # rows which have the
# same DNA accession and group
overlapRows = overlapRows[~overlapRows['DNA Accession'].isin(dupedRows['DNA Accession'])]  # Removes entries which are
# complete duplicates of one another from overlapRows, leaving only the ones that differ in everything except DNA
# Accession and group.If there are any entries in this dataframe, then searching with DNA Accession and group together
# will still result in multiple hits
protDupe = newAnn[newAnn.duplicated(subset=['Protein Accession'], keep=False)]  # find duplicate Prot. Acc.

# Cull annotations and provide the user with output detailing which entries were culled and why

# Duplicate DNA accession and group culling
cutEntries = list(overlapRows.index)
cutLines = [x+2 for x in cutEntries]  # adds 2 to index to get line # in annotation source file (+1 from counting
# from 1 instead of 0, +1 from header)
print("\033[1;31;31m The following entries (line # in the annotation file) were cut from original file because their "
      "DNA accessions and groups overlapped, making it impossible for AMR++ to read them:\033[0;31;39m")
print(cutLines)
# print("(Debug) Their index values are: ")
# print(cutEntries)
newAnn.drop(cutEntries, axis=0, inplace=True)  # removes rows with duplicate groups and DNA Accessions

# Duplicate Protein Accession culling
protDupeEntries = list(protDupe.index)
protDupeLines = [x+2 for x in protDupeEntries]
print("\033[1;31;31m The following entries (line # in the annotation file) were cut from original file because they "
      "lack or have duplicate protein accessions, making it impossible to add them to the database file:\033[0;31;39m")
print(protDupeLines)
# print("(Debug) Their index values are: ")
# print(protDupeEntries)  # only show those that were cut because of duplication, not null
newAnn.drop(protDupeEntries, axis=0, inplace=True)  # removes rows with duplicate Protein Accessions

# N/A culling
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
percentDrop = round(100 * dropTotal/len(aroIndex), 2)
print("\n Total number of entries dropped: " + str(dropTotal) + ", which is " + str(percentDrop) + "% of total entries "
                                                                                                 "\n")
#//endregion

#//region Convert multi-drug resistant class columns to "multi-drug resistant" string
# print(newAnn.loc[newAnn['class'].str.contains(';'),'class'])
newAnn.loc[newAnn['class'].str.contains(';'),'class'] = 'multi-drug resistant'

#//endregion
# TODO: Is it worth doing this if you can't also do it to the group column (because that group column is required to
#  sort annotations)? Eg. 3 entries of DNA Acc. AE004091.2 will be culled because they have the same DNA acc,
#  but different collections of multiple groups. All 3 would be found to have the same DNA acc. and groups and be culled

#//region Add "Model Name" to newAnn
#
# indexCols = ['Protein Accession', 'Model Name']  # Columns for protein accession and gene, respectively, from index file
# # Gene column is pulled from Model Name, not ARO name, because Model Name is used to build database headers
# # Create dataframe from the index file to line up each entry's gene with its gene family by protein accession
# compIndex = aroIndex[indexCols].copy()
#
# if dataframe_merge(newAnn, compIndex, doc=False, which='left_only').index.size > 0:  # Checks that there are
#     # no entries that are only in the annotation file
#     print("\033[1;31;31m ERROR: Annotation is larger than index. Ensure you have the most recent version of both \033["
#           "0;31;39m")
#     exit()
# # in the output file, "right_only" means that it is only in the index. 'left_only' should be empty, as the index
# # should have all of the entries that are in the annotation, but not vice versa. Having things in 'right_only' is
# # therefore fine, but having any in left_only means that the index does not contain all entries, and something has
# # gone wrong with your download.

# mergeCheck = dataframe_merge(newAnn, compIndex, doc=False, which='both')
# if mergeCheck[mergeCheck.duplicated(subset=['Protein Accession'], keep=False)].index.size > 0:  # check for duplicates
#     # after merge. Merge function sometimes produces duplicate entries
#     print(mergeCheck[mergeCheck.duplicated(subset=['Protein Accession'], keep=False)].index.size)
#     with pd.option_context('display.max_columns', 10):
#         print(mergeCheck[mergeCheck.duplicated(subset=['Protein Accession'], keep=False)])  # print entries
#         # duplicated by merging index and newAnn
#     print("number of duplicates:" + str(len(mergeCheck[mergeCheck.duplicated(subset=['Protein Accession'], keep='first')])))
#     print("newAnn length: " + str(len(newAnn)))
#     print("mergeCheck length: " + str(len(mergeCheck)))
#     print("ERROR: Duplicate entries detected. Exiting translator")
#     exit()
#
# newAnn = dataframe_merge(newAnn, compIndex, doc=False, which='both', ind=False)  # add 'Model Name' to newAnn,
# # merging by protein accession


#//endregion
#//endregion

#//region Database Translation

annotationAccessions = list(newAnn['DNA Accession'])  # this needs to be here because DNA Accession gets cut from
# newAnn in the next section
annotationGene = list(newAnn['Model Name'])


#//region Create AMR++-compliant header, concatenate it to newAnn, and put newAnn into the same order
typeCol = ['Drugs'] * len(newAnn.index)  # Creates a list of the string 'Drugs' with as many values as the
# annotation file has. MEGARes has a type column, but ARO (mostly) only deals with drugs
typeAnn = pd.Series(typeCol)  # Turns that list into a Dataframe so it can be concatenated to
# the new Dataframe which will contain the translated data

headerCol = newAnn['DNA Accession'].map(str) + "|" + typeAnn.map(str) + "|" + newAnn['class'].map(str) + "|" + \
            newAnn['mechanism'].map(str) + "|" + newAnn['group'].map(str) + "|" + "RequiresSNPConfirmation"
# concatenate columns to make header and add "RequiresSNPConfirmation" flag to force AMR++ to use RGI's perfect
# algorithm, because all the entries I translated are from the protein homolog model only

newAnn = pd.concat([headerCol, typeAnn, newAnn['class'], newAnn['mechanism'], newAnn['group'],
                    newAnn['Protein Accession'], newAnn['Model Name'], newAnn['DNA Accession']], axis=1)  #
# concatenates all columns that must be in the final annotation

newAnn.columns = ['header', 'type', 'class', 'mechanism', 'group', 'Protein Accession', 'Model Name', 'DNA Accession']
# rename columns to match with those of AMR++. Protein Accession and Model Name will be cut before export

#//endregion

#//region dbGenes and dbAccession Processor
dbAccessions = []
dbGenes = []
p = re.compile('(gb\|)')  # find the string 'gb|' in each database entry's header
for i in range(0,len(aroDB)):  # Create lists of DNA Accessions and genes from database headers
    dbAccessions.append(p.sub('', aroDB[i].id))  # remove 'gb|' from each entry's header
    dbAccessions[i] = dbAccessions[i][:dbAccessions[i].index('|')]  # Cuts everything after the first |, leaving only
    # the DNA accession
    dbGenes.append(aroDB[i].description)
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
dbToCull = []
overlapMessage = 'overlap culled'
protDupeMessage = 'protein accession duplication culled'
noAnnotationMessage = 'lack annotation culled'

for i in range(0, len(annotationAccessions)):
    for n in range(0, len(dbAccessions)):  # Sorts new headers into same order as database
        if annotationAccessions[i] == dbAccessions[n] and annotationGene[i] == dbGenes[n]:  # match by accession and
            # gene
            match.append([i, n])  # give list of indices of accessions that match each other
            newHeaders[n] = newAnn['header'].loc[i]  # set list to contain all translated headers in the correct
            # order for the translated database
for i in range(0, len(dbAccessions)):  # Gets indices of database entries whose annotations have been culled so that
    # those database entries can be removed later
    if dbGenes[i] in list(overlapRows['Model Name']) and dbAccessions[i] in list(overlapRows['DNA Accession']):
        # DB entries whose annotations were culled for overlapping DNA accession and gene family
        newHeaders[i] = overlapMessage
        dbToCull.append(i)
    elif dbGenes[i] in list(protDupe['Model Name']) and dbAccessions[i] in list(protDupe['DNA Accession']):
        # DB entries whose annotations were culled for having a duplicate protein accession
        newHeaders[i] = protDupeMessage
        dbToCull.append(i)
    # elif dbGenes[i] in noAnnotationGene and dbAccessions[i] in noAnnotationAccession:
    #     # DB entries with no annotation entry
    #     newHeaders[i] = noAnnotationMessage
    #     dbToCull.append(i)

# x = 1 # test value that determines which annotation/databse entry pair to print
# print(aroDB[x].description)  # print original id. If print(aroDB[x].id) matches its DNA Accession with
# # print(newHeaders[x]), the translator is creating the list of translated headers in the same order as the database,
# # allowing one-to-one indexing (aroDB[1] should have the same DNA accession as newHeaders[1]. newHeaders will have a
# # gene family where aroDB has a gene)
# print(newHeaders[x])
# # print(aroDB[x])
# print("Matching Accessions: ")
# print(match)

# Error checking before proceeding to the file writing stage
noValue = 0
# list of db entries that have not been assigned a header and were not culled, suggesting that an error has occurred
# somewhere. ProtDupe and DupedRows always have NaN for their model name, so do not need to be culled from database
errorPresent = False

for i in range(0, len(newHeaders)):  # Checks for database entries have not been assigned a header, either correctly
    # (due to culling) or erroneously
    if newHeaders[i] == overlapMessage:
        print("Database entry on line " + str((i + 1) * 2 - 1) + " Has been culled because its annotation's DNA "
                                                               "Accession and gene family overlapped with another "
                                                                 "annotation")
        # print(aroDB[i].description)
    if newHeaders[i] == noAnnotationMessage:
        print("Database entry on line " + str((i + 1) * 2 - 1) + " Has been culled because it has no corresponding annotation")
        # print(aroDB[i].description)
    if newHeaders[i] == protDupeMessage:  # This message shouldn't appear for current CARD data (Feb 2020), but will be
        # left in in case new data is added
        print("Database entry on line " + str((i + 1) * 2 - 1) + " Has been culled because its protein accession was "
                                                         "identical to another annotation")
        # print(aroDB[i].description)
    elif newHeaders[i] == 'error':  # Indicates database entries which will not be assigned a header,
        # but whose annotations were not culled, suggesting that an error in header assignment has occurred
        print("ERROR: Database entry on line " + str((i + 1) * 2 - 1) + " has not been given a value")
        print("DEBUG: index = " + str(i))
        print(aroDB[i].description)
        errorPresent = True
        noValue += 1
print("Number of unmatched entries: " + str(noValue))
print("Database entries to cull: " + str(len(dbToCull)))
print("Line numbers of culled database entries: ")
culledDB = [((line+1)*2-1) for line in dbToCull]
print(culledDB)

if errorPresent:
    print("ERROR: Some database entries are not being assigned headers, but are also not being culled. No files "
          "generated")
    exit()
#//endregion

#//region Write final files
print("Writing Annotation file")
finalCols = ['header', 'type', 'class', 'mechanism', 'group']
finalAnn = newAnn[finalCols].copy()  # drop all columns that are unneeded for
# annotation file
# with pd.option_context('display.max_columns', 5):
#     print(finalAnn)

# Write annotation file
today = date.today()  # get current date
filename = ("CARD_to_AMRplusplus_Annotation_" + today.strftime("%Y_%b_%d") + ".csv")
# Exports AMR++-ready annotation file and names it based on the present date

pd.DataFrame.to_csv(finalAnn, filename, index=False)  # exports converted annotation file as csv
print("DONE")
# Write Database file and cull problem headers
print("Writing Database file")
headerToCull = []
for i in dbToCull:
    headerToCull.append(aroDB[i].id)  # Original database headers of DB entries to be culled
    del newHeaders[i]  # Remove culled entries from newHeaders

translatedFilename = ("./CARD_to_AMRplusplus_Database_" + today.strftime("%Y_%b_%d") + ".fasta")

# Write all entries to database file that are not on the cull list
keepSeq = {}  # Hash table to contain sequences that have been given an appropriate header
for seq_record in newAroDB:  # pulls each Seq_record (header + sequence) from the source newAroDB
    header = str(seq_record.id)  # extracts the header of the seq_record as a string
    if header not in headerToCull and header not in keepSeq:  # check that the entry does not need to be culled
        keepSeq[header] = str(seq_record.seq)  # store sequence of unculled entry
    i = 0
    with open(translatedFilename, 'w') as translated:
        for header in keepSeq:
            translated.write(">" + newHeaders[i] + "\n" + keepSeq[header] + "\n")
            i += 1
print("DONE")
# //endregion
# THE END
