# Translates CARD Annotation file into one that AMR++ can read


# TODO: REMOVE Accessions, replace with unique number (followed by the accession, for readability)
# TODO: Add another round of accession culling to remove duplicates created during the newAnn creation process
# TODO: Remove Date stamping when translator is confirmed working
# TODO: Remove "multi-drug resistance" if multi-group entries cannot be removed or changed to "multi-group" somehow
# TODO: Create README file for program
# TODO: Remove private annotations from translated. Compare private list to translated list by DNA Accession and remove
#  those that match

import os
import pandas as pd
from datetime import date
from Bio import SeqIO
import re
import argparse
import csv

pd.options.display.width = 0

print("Definitions: \n"
      "Entry: Information corresponding to a single determinant contained on one line of the annotation file or two "
      "lines of the database file. \n"
      "Annotation file: Contains information used to construct new database entry headers\n"
      "Database file: Contains sequence information for each determinant \n"
      "Line number: The line number of an entry in the *original file*, not the final file. In the database file, "
      "this points to the header's line number, not the sequence's line number. \n")

# Define default filenames
aroIndexFile = 'aro_index.tsv'
aroDBFile = "nucleotide_fasta_protein_homolog_model.fasta"

#//region Parse commandline arguments to allow the user to input different filenames
parser = argparse.ArgumentParser()
parser.add_argument("-i", help='Allows user to define the filepath of the ARO index file')
parser.add_argument("-d", help='Allows user to define the filepath of the ARO database file')

args = parser.parse_args()
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


def entry_to_line(x): # Converts from index # to line #
    return (x+1)*2 - 1


# def dataframe_merge(df1, df2, doc=False, which=None, on='Protein Accession', ind=True, byIndex=False):  # Compares 2
#     # dataframes for their contents and outputs the results. set doc to true to output a csv file containing the
#     # merged dataframe. pass which='both'to check those that are in both one and the other.
#     comparison_df = df1.merge(df2,
#                               indicator=True,
#                               how='outer',
#                               on=on,
#                               left_index=byIndex,
#                               right_index=byIndex
#                               )
#     duplicates = comparison_df[comparison_df.duplicated(keep='first')].index  # finds duplicate entries created by
#     # the merge
#     comparison_df.drop(duplicates, axis=0, inplace=True)  # drops duplicate entries created by the merge
#     comparison_df.reset_index(drop=True, inplace=True)  # Reset index after dropping entries to prevent empty rows
#     if ind:
#         if which is None:
#             merged_df = comparison_df[comparison_df['_merge'] != 'both']  # returns only the entries which differ
#         else:
#             merged_df = comparison_df[comparison_df['_merge'] == which]
#     else:
#         merged_df = comparison_df[comparison_df['_merge'] == 'both']
#         # merged_df.drop(['_merge'], axis=1, inplace=True)
#     if doc:
#         merged_df.to_csv('diff.csv')
#     return merged_df


#//endregion

# Import files

# privateAnn = pd.read_csv('private_models.csv')  # Read in annotations that are not used by CARD in order to remove
# them from the translated annotation list
aroIndex = pd.read_csv(aroIndexFile, sep='\t')  # Read index file in order to generate annotations
aroDB = list(SeqIO.parse(aroDBFile, "fasta"))  # Read in CARD Database file as a list of seqRecord objects
newAroDB = SeqIO.parse(aroDBFile, "fasta")  # Read in CARD database as a Seqrecord generator, instead of a list,
# so that it can be overwritten later

#//region Annotation Translation

# Create new DataFrame containing AMR++-relevant data
aroCols = ['DNA Accession', 'Drug Class', 'Resistance Mechanism', 'AMR Gene Family',
           'Model Name']  # Important
# columns from ARO data.
newAnn = aroIndex[aroCols].copy()  # Creates a Dataframe containing the Columns from the ARO annotation
# file. ARO's columns (branches) are in a different order than MEGARes, so this reorders them.
# Also Add DNA accession to allow the database to be searched for matching entries - fills the same spot as Meg_###

newAnn.columns = ['DNA Accession', 'class', 'mechanism', 'group', 'Model Name']  # sets names of columns of new
# annotation file

# Replace spaces with underscores because AMR++ can't parse spaces
newAnn['group'] = newAnn['group'].str.replace(' ', '_')
newAnn['mechanism'] = newAnn['mechanism'].str.replace(' ', '_')
newAnn['class'] = newAnn['class'].str.replace(' ', '_')
newAnn['Model Name'] = newAnn['Model Name'].str.replace(' ', '_')

for i in range(0,len(aroDB)):
    aroDB[i].description = aroDB[i].description.replace(' ', '_')

#//region -PRIME notation conversion
# Removes parentheses in gene family and converts ' to -PRIME and '' to -DPRIME. eg. AAC(6') = AAC6-PRIME
newAnn['group'] = newAnn['group'].str.replace('(?<! )\(', '', regex=True)
# replaces open paren only when within a word
newAnn['group'] = newAnn['group'].str.replace('\'\'\)$', '-DPRIME', regex=True)
newAnn['group'] = newAnn['group'].str.replace('\'\)$', '-PRIME', regex=True)  #
# Replace ' or '' when next to a close paren and at the end of a line (single gene family)
newAnn['group'] = newAnn['group'].str.replace('\'\'\);', '-DPRIME;', regex=True)
newAnn['group'] = newAnn['group'].str.replace('\'\);', '-PRIME;', regex=True)
# Replace ' or '' when next to a close paren and next to a semicolon (multiple gene families)
newAnn['group'] = newAnn['group'].str.replace('(?<=\d)\)$', '', regex=True)
newAnn['group'] = newAnn['group'].str.replace('(?<=\d)\);', ';', regex=True)
# replaces ) when preceded by a number and (at the end of a string or semicolon-separated)
#//endregion


#//region Create AMR++-compliant header, concatenate it to newAnn, and put newAnn into the same order
typeCol = ['Drugs'] * len(newAnn.index)  # Creates a list of the string 'Drugs' with as many values as the
# annotation file has. MEGARes has a type column, but ARO (mostly) only deals with drugs
typeAnn = pd.Series(typeCol)  # Turns that list into a Dataframe so it can be concatenated to
# the new Dataframe which will contain the translated data

headerCol = newAnn['DNA Accession'].map(str) + "|" + typeAnn.map(str) + "|" + newAnn['class'].map(str) + "|" + \
            newAnn['mechanism'].map(str) + "|" + newAnn['group'].map(str)
# concatenate columns to make header

# headerCol = newAnn['DNA Accession'].map(str) + "|" + typeAnn.map(str) + "|" + newAnn['class'].map(str) + "|" + \
#             newAnn['mechanism'].map(str) + "|" + newAnn['group'].map(str) + "|" + "RequiresSNPConfirmation"
# RequiresSNPConfirmation is only required for Protein Variant Model, not PHM. Swap this out if adding PVM to
# translator.

newAnn = pd.concat([headerCol, typeAnn, newAnn['class'], newAnn['mechanism'], newAnn['group'], newAnn['Model Name'],
                    newAnn['DNA Accession']], axis=1)  #
# concatenates all columns that must be in the final annotation

newAnn.columns = ['header', 'type', 'class', 'mechanism', 'group', 'Model Name', 'DNA Accession']
# rename columns to match with those of AMR++. Protein Accession and Model Name will be cut before export

#//endregion


#//region Unsearchable annotation checking and culling

# Check for annotations that cannot be searched
dupedRows = newAnn[newAnn.duplicated(subset=['header'],
                                     keep=False)].copy()  # Only want to remove those that are duplicates in the
# areas that were used to create the header as that is what will be attached to the database

# Rows that are just duplicate annotations
overlapRows = newAnn[newAnn.duplicated(subset=['DNA Accession', 'Model Name'],
                                       keep=False)].copy()  # rows which have the same DNA accession and gene

overlapRows = overlapRows[~overlapRows['DNA Accession'].isin(dupedRows['DNA Accession'])]  # Removes entries which are
# complete duplicates of one another from overlapRows, leaving only the ones that differ in everything except DNA
# Accession and gene. If there are any entries in this dataframe, then searching with DNA Accession and gene together
# will still result in multiple hits

#
# protDupe = newAnn[newAnn.duplicated(subset=['Protein Accession'],
#                                     keep=False)].copy()  # find duplicate Prot. Acc.

# Cull annotations and provide the user with output detailing which entries were culled and why

# Duplicate DNA accession and group culling
overlapEntries = list(overlapRows.index)
overlapLines = [x+2 for x in overlapEntries]  # adds 2 to index to get line # in annotation source file (+1 from counting
# from 1 instead of 0, +1 from header)
print("\033[1;31;31m The following entries (line # in the annotation file) were cut from original file because their "
      "DNA accessions and genes overlapped, causing them to conflict when being attached to a database entry:\033["
      "0;31;39m")
print(overlapLines)
# print("(Debug) Their index values are: ")
# print(overlapEntries)
newAnn.drop(overlapEntries, axis=0, inplace=True)  # removes rows with duplicate groups and DNA Accessions

# Duplicate Protein Accession culling
# protDupeEntries = list(protDupe.index)
# protDupeLines = [x+2 for x in protDupeEntries]
# print("\033[1;31;31m The following entries (line # in the annotation file) were cut from original file because they "
#       "lack or have duplicate protein accessions, making it impossible to add them to the database file:\033[0;31;39m")
# print(protDupeLines)
# # print("(Debug) Their index values are: ")
# # print(protDupeEntries)  # only show those that were cut because of duplication, not null
# newAnn.drop(protDupeEntries, axis=0, inplace=True)  # removes rows with duplicate Protein Accessions

# N/A culling
nullEntries = newAnn[newAnn.isna().any(axis=1)].index.tolist()  # finds indices of any entry with any n/a cell
nullLines = [x+2 for x in nullEntries]
nullGene = newAnn["Model Name"].loc[nullEntries].tolist() # List of genes of entries with null values
nullAccession = newAnn["DNA Accession"].loc[nullEntries].tolist() # List of DNA Accessions of entries with null
# values. These two together will be used to search the database for entries that need to be culled

print("\033[1;31;31m The following entries (line # in the annotation file) were cut from original file because they "
      "lack a DNA Accession, drug class, mechanism, or gene family, making it impossible to add "
      "them to the database file:\033[0;31;39m")
print(nullLines)
# print("(Debug) Their index values are: ")
# print(Diff(nullEntries, protDupeEntries))
newAnn.dropna(how='any', axis=0, inplace=True)  # removes rows with n/a in any field
newAnn.reset_index(drop=True, inplace=True)  # Reset index after dropping entries to prevent empty rows

dropTotal = len(nullEntries)
percentDrop = round(100 * dropTotal/len(aroIndex), 2)
print("\n Total number of annotation entries dropped: " + str(dropTotal) + ", which is " + str(percentDrop) + "% of "
                                                                                                            "total entries "
                                                                                                 "\n")
#//endregion
#//endregion

#//region Database Translation

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
    dbGenes[i] = splitGroup[5]  # adds the 'gene' section of the database header into the dbGenes list
    dbGenes[i] = dbGenes[i][:dbGenes[i].index('_[')]  # cuts species information, leaving only the gene
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


#//region Prepare and match new headers to overwrite database headers
match = []
newHeaders = ['error'] * len(dbAccessions)  # Create list that will contain all translated headers in the correct order
# for the translated database. If one is not filled in, it will be listed as "error". Must create it at full length
# to allow it to be filled out of order, as the database and annotation files are in different orders.
dbToCull = []
granularityMessage = 'granularity culled'
overlapMessage = 'overlap culled'
protDupeMessage = 'protein accession duplication culled'
noAnnotationMessage = 'lack annotation culled'
nullEntryMessage = 'null entry culled'
wrongHeaders = [protDupeMessage, nullEntryMessage, noAnnotationMessage, overlapMessage, granularityMessage]

annotationAccessions = list(newAnn['DNA Accession'])
annotationGene = list(newAnn['Model Name'])

# Assign new headers
for i in range(0, len(annotationAccessions)):
    for n in range(0, len(dbAccessions)):  # Sorts new headers into same order as database
        # if newAnn['header'].loc[i] in newHeaders:
        #     newHeaders[n] = granularityMessage
        #     break  # Break to prevent it from overwriting the entire newHeader with granularityMessage
        if annotationAccessions[i] == dbAccessions[n] and annotationGene[i] == dbGenes[n]: # match by accession
            # and gene
            match.append([i, n])  # give list of indices of accessions that match each other. DEBUG only.
            newHeaders[n] = newAnn['header'].loc[i]  # set list to contain all translated headers in the correct
            # order for the translated database

granularityCulled = []
for i in range(len(newHeaders) - 1, -1, -1):  # go through newHeaders in reverse order and remove duplicate headers
    if newHeaders.count(newHeaders[i]) > 1 and newHeaders[i] != granularityMessage:
        granularityCulled.append(newHeaders[i])
        newHeaders[i] = granularityMessage

# Get indices of database entries whose annotations have been culled so that those database entries can be removed
for i in range(0, len(dbAccessions)):

    if dbGenes[i] in list(overlapRows['Model Name']) and dbAccessions[i] in list(overlapRows['DNA Accession']):
        # DB entries whose annotations were culled for overlapping DNA accession and gene
        newHeaders[i] = overlapMessage
        dbToCull.append(i)

    elif dbGenes[i] in nullGene and dbAccessions[i] in nullAccession:
        newHeaders[i] = nullEntryMessage
        dbToCull.append(i)

    elif newHeaders[i] == granularityMessage:
        dbToCull.append(i)
#//endregion

#//region Error checking before proceeding to the file writing stage
noValue = 0  # list of db entries that have not been assigned a header and were not culled, suggesting that an error
# has occurred somewhere. ProtDupe and DupedRows always have NaN for their model name, so do not need to be culled
# from database
unassignedHeaders = []

#TODO: Instead of printing each one indvidually and outputting an illegible mess, just output a list of all the
# culled entries in each category (culled for overlap, no annotation, protdupe, granularityCulled, etc)

for i in range(0, len(newHeaders)):  # Checks for database entries have not been assigned a header, either correctly
    # (due to culling) or erroneously
    if newHeaders[i] == overlapMessage:
        print("Database entry on line " + str(entry_to_line(i)) + " Has been culled because its annotation's DNA "
                                                               "Accession and gene family overlapped with another "
                                                                 "annotation")
        # print(aroDB[i].description)
    elif newHeaders[i] == noAnnotationMessage:
        print("Database entry on line " + str(entry_to_line(i)) + " Has been culled because it has no corresponding annotation")
        # print(aroDB[i].description)
    elif newHeaders[i] == protDupeMessage:  # This message shouldn't appear for current CARD data (Feb 2020),
        # but will be left in in case new data is added
        print("Database entry on line " + str(entry_to_line(i)) + " Has been culled because its protein accession was "
                                                         "identical to another annotation")
        # print(aroDB[i].description)
    elif newHeaders[i] == nullEntryMessage:
        print("Database entry on line " + str(entry_to_line(i)) + " Has been culled because its annotation contained a "
                                                                  "null value")
    elif newHeaders[i] == granularityMessage:
        print("Database entry on line " + str(entry_to_line(i)) + " Has been culled because its annotation was "
                                                                  "already assigned to another sequence. AMR++ only "
                                                                  "allows one header per sequence, and the conversion "
                                                                  "from gene to gene families results in a loss of "
                                                                  "granularity.")
        # print(aroDB[i].description)
    elif newHeaders[i] == 'error':  # Indicates database entries which will not be assigned a header,
        # but whose annotations were not culled, suggesting that an error in header assignment has occurred
        print("ERROR: Database entry on line " + str(entry_to_line(i)) + " has not been given a value")
        print("DEBUG: index = " + str(i))
        print(aroDB[i].description)
        unassignedHeaders.append(aroDB[i].description)
        noValue += 1
print("Number of unmatched entries: " + str(noValue))
percentCulled = round(len(dbToCull)/len(aroDB), 3)*100
print("Database entries to cull: " + str(len(dbToCull)) + ", or " + str(percentCulled) + "% of the CARD "
                                                                                         "Homolog Model Database")
print("Line numbers of culled database entries (In original database file): ")
culledDB = [entry_to_line(line) for line in dbToCull]
print(culledDB)

if len(unassignedHeaders) > 0:
    print("ERROR: Some database entries are not being assigned headers, but are also not being culled. No files "
          "generated. They are:")
    print(unassignedHeaders)
    exit()

# Check that there is only one of each kind of newHeader. Every database header must be unique.
DuplicatePrinted = False
for i in newHeaders:
    if newHeaders.count(i) > 1 and i not in wrongHeaders:
        if not DuplicatePrinted:
            print("ERROR: Multiple copies of the same translated header are being created. No files have been "
                  "generated. The following headers are being duplicated: ")
            DuplicatePrinted = True
        print(i)
if DuplicatePrinted:
    exit()

#//endregion
#//endregion

#//region Write final files

# Write annotation file. Exports AMR++-ready annotation file and names it based on the present date

translatedPath = './translations/'
today = date.today()  # get current date
if not os.path.exists(translatedPath):  # Creates a folder for translted files to go into if one is not already present
    os.mkdir(translatedPath)
annotationFilename = (translatedPath + "CARD_to_AMRplusplus_Annotation_" + today.strftime("%Y_%b_%d") + ".csv")
translatedFilename = (translatedPath + "CARD_to_AMRplusplus_Database_" + today.strftime("%Y_%b_%d") + ".fasta")
culledDBFilename = 'Overlap_culled_DB.csv'

print("Writing Annotation file")

finalCols = ['header', 'type', 'class', 'mechanism', 'group']
finalAnn = newAnn[finalCols].copy()  # drop all columns that are unneeded for
# annotation file

# Export converted annotation file as csv
pd.DataFrame.to_csv(finalAnn, annotationFilename, index=False)
print("DONE")

# Write Database file and cull problem headers
print("Writing Database file")
headerToCull = []

for i in sorted(dbToCull, reverse=True):  # Deletes entries in reverse order to prevent changing the indices of the
    # list as you delete them
    headerToCull.append(aroDB[i].id)  # Original database headers of DB entries to be culled
    del newHeaders[i]  # Remove culled entries from newHeaders

if any(item in wrongHeaders for item in newHeaders):
    print("ERROR: culled headers are not being dropped from newHeaders and they are being sent to the database. "
          "Annotation file already generated, but no database file generated.")
    exit()


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

print("Writing Translated versions of database entries culled due to annotation overlap (reducing granularity) into "
      "csv... ") # Print
print("Filename: " + culledDBFilename)
# duplicate headers to allow comparison of their families to MEGARes
with open(culledDBFilename, 'w', newline='') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(granularityCulled)

print("DONE")

# //endregion

# THE END
