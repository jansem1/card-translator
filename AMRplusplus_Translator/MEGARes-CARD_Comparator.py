# - Do -PRIME notation conversion on CARD data
# - Connect each sequence to its CARD/MEGARes gene family - 2 dataframes
# Merge the dataframes along the sequence to produce a final dataframe:
# Dataframe: [Sequence, CARD gene fam, MEGARes group]
# - Figure out where the families and groups overlap - bar charts

# TODO: Switch to translated CARD database? That's what's actually being compared

import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

#//region Define functions


def importFasta(file):  # Pandas doesn't have its own fasta importer, so we bring in fasta data with SeqIO and
    # convert it to a dataframe
    with open(file) as fasta_file:  # Will close handle cleanly
        identifiers = []
        sequences = []
        for title, sequence in SimpleFastaParser(fasta_file):
            identifiers.append(title)  # First word is ID
            sequences.append(sequence)
    d = list(zip(identifiers,sequences))
    return pd.DataFrame(d, columns=['header', 'sequence'])


def get_groups (data, grouploc, species=False):  # Extracts group from header
    groups = data['header'].str.split('|')
    f = lambda x: x[grouploc]  # Megares and CARD have their groups at different locations in the
    # header, so the group's location must be specified.
    groups = groups.apply(f)
    groups.columns = ['group']

    if species:  # Only CARD has species info
        g = lambda x: x[:x.index(' [')]  # removes species name. Can't just remove by space because some group names
        # have spaces in them
        groups = groups.apply(g)

    out = data.merge(groups, left_index=True, right_index=True)
    out.columns = ['header', 'sequence', 'group']  # The merge changes the header and group column names for some
    # reason, so they have to be changed back
    return out


def prime_notation(data):  # Converts from (#') to #-PRIME
    # Removes parentheses in gene family and converts ' to -PRIME and '' to -DPRIME. eg. AAC(6') = AAC6-PRIME
    data['group'] = data['group'].str.replace('(?<! )\(', '')
    # replaces open paren only when within a word
    data['group'] = data['group'].str.replace('\'\'\)$', '-DPRIME')
    data['group'] = data['group'].str.replace('\'\)$', '-PRIME')  #
    # Replace ' or '' when next to a close paren and at the end of a line (single gene family)
    data['group'] = data['group'].str.replace('\'\'\);', '-DPRIME;')
    data['group'] = data['group'].str.replace('\'\);', '-PRIME;')
    # Replace ' or '' when next to a close paren and next to a semicolon (multiple gene families)
    data['group'] = data['group'].str.replace('(?<=\d)\)$', '')
    data['group'] = data['group'].str.replace('(?<=\d)\);', ';')
    # replaces ) when preceded by a number and (at the end of a string or semicolon-separated)
    return data

#//endregion


megDataFile = 'megares_full_database_v2.00.fasta'
cardDataFile = 'nucleotide_fasta_protein_homolog_model.fasta'

megData = importFasta(megDataFile)
cardData = importFasta(cardDataFile)
# Columns are now 0: Header, 1: Sequence

# Pull out group/family and append it to the end of the dataframe
cardData = get_groups(cardData, 5, species=True)
megData = get_groups(megData, 4)
# Columns are now 0: Header, 1: Sequence, 2: group

cardData = prime_notation(cardData)
pd.set_option('display.max_rows', None, 'display.max_columns', None)
print(cardData['header'].loc[2507])