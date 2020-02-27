# - Do -PRIME notation conversion on CARD data
# - Connect each sequence to its CARD/MEGARes gene family - 2 dataframes
# Merge the dataframes along the sequence to produce a final dataframe:
# Dataframe: [Sequence, CARD gene fam, MEGARes group]
# - Figure out where the families and groups overlap - bar charts

import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from datetime import date

# Dataframe display options
pd.options.display.width = 0
pd.options.display.max_rows = 50
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


def get_bins (data, binloc, source, species=False):  # Extracts group from header

    bins = data['header'].str.split('|')
    f = lambda x: x[binloc]  # Megares and CARD have their bins at different locations in the
    # header, so the group's location must be specified.
    bins = bins.apply(f)
    bins.columns = [source + '_bins']
    if species:  # Only CARD's original database has species info
        g = lambda x: x[:x.index(' [')]  # removes species name. Can't just remove by space because some group names
        # have spaces in them
        bins = bins.apply(g)

    out = data.merge(bins, left_index=True, right_index=True)
    out.columns = ['header', 'sequence', source + '_bins']  # The merge changes the header and group column names for
    # some
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

# Pull out group/family and append it to the end of the dataframe
cardData = get_bins(cardData, 4, 'card')
megData = get_bins(megData, 4, 'meg')
# Columns are now 0: Header, 1: Sequence, 2: {database}_bin

x = 0
for i in cardData.index:
    if cardData['sequence'].loc[i] in megData['sequence'].tolist():
        x+=1

cardInMeg = (str(round(100 * x/len(cardData.index), 2)))
megInCard = (str(round(100 * x/len(megData.index), 2)))
print("Percent of CARD homolog model sequences in MEGARes: " + cardInMeg + "%")
print("Percent of MEGARes sequences in CARD homolog model: " + megInCard + "%")

#//region Find all group/family matches

mergedDatabases = pd.merge(left=cardData, right=megData, on='sequence')
mergedDatabases.rename(columns={'header_x': 'card_header', 'header_y': 'meg_header'}, inplace=True)

#//endregion

print(mergedDatabases)

# TODO: Create 2 dataframes - 1 contains one instance of each MEG group and all the families that fall into it,
# 2 contains one instance of each CARD family and all the groups that fall into it.
# TODO: Delete duplicate bins (x) and duplicate bins of bins (y)

# for  i in mergedDatabases.index:
goneThrough = []
card_binsOfBins = []
for i in mergedDatabases.index:
    if mergedDatabases['card_bins'].loc[i] not in goneThrough:
        print(mergedDatabases['card_bins'].loc[i])
        print(mergedDatabases['meg_bins'].loc[mergedDatabases['card_bins'] == mergedDatabases['card_bins'].loc[i]].tolist())
    goneThrough.append(mergedDatabases['card_bins'].loc[i])