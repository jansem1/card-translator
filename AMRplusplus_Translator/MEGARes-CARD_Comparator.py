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

def bin_overlap(data, left, right):
    # TODO: THIS IS BROKEN. SPITS OUT WRONG MATCHES SOMTIMES, BUT NOT ALWAYS Eg. HERA is not included
    left_goneThrough = []
    # right_goneThrough = []
    left_binsOfBins = []
    right_containedBins = []
    # rightdoubles = []
    # leftdoubles = []
    for i in data.index:
        leftFam = data[left + '_bins'].loc[i]
        rightGroup = data[right + '_bins'].loc[data[left + '_bins'] == leftFam].tolist()  # Finds
        # all
        # instances of MEGARes groups that overlap with each CARD family
        if leftFam not in left_goneThrough:
            left_binsOfBins.append(leftFam)
            right_containedBins.append(rightGroup)
            # if rightGroup in right_goneThrough:
            #     rightdoubles.append(rightGroup)
            #     leftdoubles.append(leftFam)
            # if rightGroup not in right_goneThrough:
            #     right_goneThrough.append(rightGroup)
        left_goneThrough.append(data[left + '_bins'].loc[i])
    df = {left: left_binsOfBins, right: right_containedBins}
    return pd.DataFrame(df)


removeDuplicates = lambda x: list(dict.fromkeys(x))
# TODO: This breaks apart card families into unique characters. Stop that.

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

# TODO: Create 2 dataframes - 1 contains one instance of each MEG group and all the families that fall into it,
#  2 contains one instance of each CARD family and all the groups that fall into it.
# TODO: Delete duplicate bins (x) and duplicate bins of bins (y)

# print(mergedDatabases.loc[mergedDatabases['meg_bins'] == 'PNGM'])
print(type(mergedDatabases['meg_bins'] == 'PNGM'))
# exit()

cardOverlap = bin_overlap(mergedDatabases, 'card', 'meg')
cardOverlap['meg'] = cardOverlap['meg'].apply(removeDuplicates)
print(cardOverlap)
# exit()
# print(type(cardOverlap['meg'].loc[216]))
# print(cardOverlap['meg'].loc[216])
# print(cardOverlap['meg'].loc[216] == ['PNGM'])
# # exit()

# print(type(cardOverlap['meg'] == 'PNGM'))
# print(cardOverlap.loc[cardOverlap['meg'] == ['PNGM']])
# exit()
# TODO: bin_overlap is causing issues and unclear if it's aligning properly. Try to check by looking for all the
#  entries in cardOverlap that contain "HERA" in the 'meg' column. That works fine in mergedDatabases, but not in
#  cardOverlap.
megOverlap = bin_overlap(mergedDatabases, 'meg', 'card')
megOverlap['card'] = megOverlap['card'].apply(removeDuplicates)
print(megOverlap)
exit()
# print(len(megOverlap['card']))
print(megOverlap.loc[megOverlap['card'].map(len) > 1])  # Finds MEG entries with more than one family
