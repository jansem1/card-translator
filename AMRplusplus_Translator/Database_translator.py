# Translates CARD database's header into a header that AMR++ can read. MUST BE RUN AFTER ANNOTATION TRANSLATOR TO GET
# HEADER DATA

# TODO: Compare DNA Accession from header to DNA Accession from database file, take Annotation header and replace
#  matching DB headers
# TODO #  3. If the DNA Accessions match, replace the aroDB string with its annotation equivalent (STARTING WITH A '>'!)
#  4. Write new FASTA file
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
# filepath = easygui.fileopenbox('test') # Alternative to current-day file system. Currently broken by issues
# with tkinter
today = date.today()
filename = ("CARD_to_AMRplusplus_Annotation_" + today.strftime("%Y_%b_%d") + ".csv")  # Opens annotation file that was created today
annotationSource = pd.read_csv(filename)   # Reads file
# Extract header from translated annotation instead of original index file in order to minimize dependencies

# print(annotationSource['header'].loc[1][:annotationSource['header'].loc[1].index('|')])
annotationAccessions = []
for i in range(0,len(annotationSource)):
    annotationAccessions.append(annotationSource['header'].loc[i][:annotationSource['header'].loc[i].index('|')])
    # print(annotationAccessions[i])

dbAccessions = []
p = re.compile('(gb\|)')  # find the string 'gb|'
for i in range(0,len(aroDB)):
    dbAccessions.append(p.sub('', aroDB[i].id))  # remove 'gb|' from each entry's header
    # aroDB[i].id = p.sub('', aroDB[i].id)
    dbAccessions[i] = dbAccessions[i][:dbAccessions[i].index('|')]  # Cuts everything after the first |, leaving only the DNA accession
    # dbAccessions.append(aroDB[i].id[:aroDB[i].id.index('|')])
    # print(dbAccessions[i])

x = 0 # test value that determines which annotation/databse entry pair to print
print(aroDB[x].id) # print original id. If print(aroDB[n].id) matches with print(newAroDB[n]), the translator is working
# creating the list of translated headers in the right order

match = []
newAroDB = ['void'] * len(dbAccessions)  # Create list that will contain all translated headers in the correct order for
# the translated database
for i in range(0, len(annotationAccessions)):
    for n in range(0, len(dbAccessions)):
        if annotationAccessions[i] == dbAccessions[n]:
            match.append([i, n])  # give list of accessions that match each other
            newAroDB[n] = annotationSource['header'].loc[i]  # set list to contain all translated headers in the correct
            # order for the translated database
# TODO: Write new FASTA with headers replaced with newAroDB
print(newAroDB[x])
# print(aroDB[0])
# print("Matching Accessions: ")
# print(match)
