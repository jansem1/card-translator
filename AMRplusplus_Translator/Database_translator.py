# Translates CARD database's header into a header that AMR++ can read. MUST BE RUN AFTER ANNOTATION TRANSLATOR TO GET
# HEADER DATA

# TODO: import CARD database (nucleotide_fasta_...), extract the header, replace with new headers. Will need to find
#  each header via DNA accession and replace it
# TODO: Compare DNA Accession from header to DNA Accession from database file, take Annotation header and replace
#  matching DB headers
# TODO go through all AroDB items, remove characters until you've removed the first | (they all start with >gb|)
#  2. Compare the AroDB string, up to the next |, to the annotation list (up to its first |). Use "in" operator and
#  "map", alongside re.sub().
#  3. If they match, replace the AroDB string with its annotation equivalent (STARTING WITH A '>'!)
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

AroDB = list(SeqIO.parse("nucleotide_fasta_protein_homolog_model.fasta", "fasta"))
# headerPath = easygui.fileopenbox('test')
today = date.today()
headerSource = pd.read_csv('aro_categories_index.tsv', sep='\t')  # Opens
# annotation file that was created today
print(headerSource['header'])
p = re.compile ('(gb\|)')  # find the string 'gb|'
for i in range(0,len(AroDB)):
    # print(AroDB[i].id)
    AroDB[i].id = p.sub('', AroDB[i].id)  # remove 'gb|' from each entry's header
    # print(AroDB[i].id)
    AroDB[i].id = AroDB[i].id[:AroDB[i].id.index('|')]  # Cuts everything after the first |, leaving only the DNA accession
    # print(AroDB[i].id)

for i in range(0,len(headerSource['header']))
    headerSource[i]

# prePipe = re.compile('(.+?(?=\|)')  # find everything before the first pipe


