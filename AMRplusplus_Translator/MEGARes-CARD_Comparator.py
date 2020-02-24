# - Add gene families/groups to dataframe
# - Do -PRIME notation conversion on CARD data
# - Connect each sequence to its CARD/MEGARes gene family - 2 dataframes
# Merge the dataframes along the sequence to produce a final dataframe:
# Dataframe: [Sequence, CARD gene fam, MEGARes group]
# - Figure out where the families and groups overlap - bar charts

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

def importFasta(file):  # Pandas doesn't have its own fasta importer, so we bring in fasta data with SeqIO and
    # convert it to a dataframe
    with open(file) as fasta_file:  # Will close handle cleanly
        identifiers = []
        sequences = []
        for title, sequence in SimpleFastaParser(fasta_file):
            identifiers.append(title.split(None, 1)[0])  # First word is ID
            sequences.append(sequence)
    d = list(zip(identifiers,sequences))
    return pd.DataFrame(d, columns=['header', 'sequence'])


def get_groups (data, grouploc):  # Extracts group from header
    groups = data['header'].str.split('|')
    f = lambda x: x[grouploc]  # Megares and CARD have their groups at different locations in the header,
    # so the group's location must be specified
    groups = groups.apply(f)
    groups.columns = ['group']
    return groups

def add_groups (data, groups):
    out = data.merge(groups, left_index=True, right_index=True)
    out.columns = ['header', 'sequence', 'group']  # The merge changes the header and group column names for some
    # reason, so they have to be changed back
    return out

megDataFile = 'megares_full_database_v2.00.fasta'
cardDataFile = 'nucleotide_fasta_protein_homolog_model.fasta'

megData = importFasta(megDataFile)
cardData = importFasta(cardDataFile)

# Pull out group/family
cardGroups = get_groups(cardData, 5)
megGroups = get_groups(megData, 4)

# Add that family in its own column
cardData = add_groups(cardData, cardGroups)
megData = add_groups(megData, megGroups)


print(cardData)