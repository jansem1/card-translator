# - Do -PRIME notation conversion on CARD data
# - Connect each sequence to its CARD/MEGARes gene family - 2 dataframes
# Merge the dataframes along the sequence to produce a final dataframe:
# Dataframe: [Sequence, CARD gene fam, MEGARes group]
# - Figure out where the families and groups overlap - bar charts

import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from datetime import date

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


def get_bins (data, binloc, species=False):  # Extracts group from header

    bins = data['header'].str.split('|')
    f = lambda x: x[binloc]  # Megares and CARD have their bins at different locations in the
    # header, so the group's location must be specified.
    bins = bins.apply(f)
    bins.columns = ['bins']
    if species:  # Only CARD's original database has species info
        g = lambda x: x[:x.index(' [')]  # removes species name. Can't just remove by space because some group names
        # have spaces in them
        bins = bins.apply(g)

    out = data.merge(bins, left_index=True, right_index=True)
    out.columns = ['header', 'sequence', 'bins']  # The merge changes the header and group column names for some
    # reason, so they have to be changed back
    return out


#//endregion


megDataFile = 'megares_full_database_v2.00.fasta'

translatedPath = './translations/'
today = date.today()
cardDataFile = (translatedPath + "CARD_to_AMRplusplus_Database_" + today.strftime("%Y_%b_%d") + ".fasta")

megData = importFasta(megDataFile)
cardData = importFasta(cardDataFile)
# Columns are now 0: Header, 1: Sequence

print(cardData['header'].loc[0])

# Pull out group/family and append it to the end of the dataframe
cardData = get_bins(cardData, 4)
megData = get_bins(megData, 4)
# Columns are now 0: Header, 1: Sequence, 2: bin

# Create column for corresponding group/family
cardData['in group'] = ""
megData['in family'] = ""
# Columns are now 0: Header, 1: Sequence, 2: bin, 3: in {bin}

x = 0
for i in cardData.index:
    if cardData['sequence'].loc[i] in megData['sequence'].tolist():
        x+=1

cardInMeg = (str(round(100 * x/len(cardData.index), 2)))
megInCard = (str(round(100 * x/len(megData.index), 2)))
print("Percent of CARD homolog model sequences in MEGARes: " + cardInMeg + "%")
print("Percent of MEGARes sequences in CARD homolog model: " + megInCard + "%")

#//region Find all group/family matches
# cardData['in group'] = np.where(cardData['sequence'] == megData['sequence'], 'True', 'False')
# megData['in family'] = np.where(cardData['sequence'] == megData['sequence'], 'True', 'False')

for i in megData.index:
    print(i)
    for n in cardData.index:
        if cardData['sequence'].loc[n] == megData['sequence'].loc[i]:
            cardData['in group'].loc[n] = megData['bins'].loc[i]
            megData['in family'].loc[i] = cardData['bins'].loc[n]
# TODO: Find a way to do this with generators. This is way too slow.
# //endregion

print(cardData)