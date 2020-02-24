# - Get Megares and CARD databases
# - From the headers, get their gene families/groups
# - Do -PRIME notation conversion on CARD data
# - Connect each sequence to its CARD/MEGARes gene family - 2 dataframes
# Merge the dataframes along the sequence to produce a final dataframe:
# Dataframe: [Sequence, CARD gene fam, MEGARes group]
# - Figure out where the families and groups overlap - bar charts

import pandas as pd
from Bio import SeqIO

megDataFile = 'megares_full_database_v2.00.fasta'
cardDataFile = 'nucleotide_fasta_protein_homolog_model.fasta'