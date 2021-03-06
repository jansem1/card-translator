# Translates CARD database's header into a header that AMR++ can read. MUST BE RUN AFTER ANNOTATION TRANSLATOR TO GET
# HEADER DATA

# TODO: CHANGE SEARCH SYSTEM. MULTIPLE ANNOTATIONS CAN HAVE THE SAME DNA ACCESSION
# TODO: Fix entries 1301, 2086, and 2128 not being converted properly
# TODO: Remove database entries that were removed by
import numpy as np
import pandas as pd
from Bio import SeqIO
import re
from datetime import date
# from tkinter import *
# import easygui as eg

# def headerComp (dbIn, headerSource):  # Read in and compare database and annotation headers
    # print(dbIn in headerSource)

aroDB = list(SeqIO.parse("nucleotide_fasta_protein_homolog_model.fasta", "fasta"))
# filepath = eg.fileopenbox('test') # Alternative to current-day file system. Currently broken by issues
# with tkinter
today = date.today()
filename = ("CARD_to_AMRplusplus_Annotation_" + today.strftime("%Y_%b_%d") + ".csv")  # Opens annotation file that was
# created today
annotationSource = pd.read_csv(filename)   # Reads file
# Extract header from translated annotation instead of original index file in order to minimize dependencies

# print(annotationSource['header'].loc[1][:annotationSource['header'].loc[1].index('|')])

annotationAccessions = []
annotationGroup = []
for i in range(0,len(annotationSource)):  # Create list of DNA Accession from translated annotation file
    annotationAccessions.append(annotationSource['header'].loc[i][:annotationSource['header'].loc[i].index('|')])  #
    # makes list of accessions from annotations
    annotationGroup.append(annotationSource['group'].loc[i])
    # makes list of groups from annotations
    # print(annotationAccessions[i])
# print(annotationGroup)

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
            newHeaders[n] = annotationSource['header'].loc[i]  # set list to contain all translated headers in the correct
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

#//region Write file
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
#//endregion
