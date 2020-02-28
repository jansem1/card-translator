# - Do -PRIME notation conversion on CARD data
# - Connect each sequence to its CARD/MEGARes gene family - 2 dataframes
# Merge the dataframes along the sequence to produce a final dataframe:
# Dataframe: [Sequence, CARD gene fam, MEGARes group]
# - Figure out where the families and groups overlap - bar charts

#TODO: Remove Date stamping when translator is confirmed
# TODO: Get list of families in CARD that aren't in MEGARes and vice versa

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


def get_bins (data, binloc, source, species=False):  # Extracts family/group from header
    bins = data['header'].str.split('|')
    f = lambda x: x[binloc]  # Legacy. Megares and original CARD can have their bins at different locations in the
    # header, so the group's location must be specified. Tranlsated CARD has its bins in the same location,
    # so this is not currently needed.
    bins = bins.apply(f)
    bins.columns = [source + '_bins']
    if species:  # Legacy. Only CARD's original database has species info. This is not necessary for translated CARD or
        # MEGARes
        g = lambda x: x[:x.index(' [')]  # removes species name. Can't just remove by space because some group names
        # have spaces in them
        bins = bins.apply(g)

    out = data.merge(bins, left_index=True, right_index=True)
    out.columns = ['header', 'sequence', source + '_bins']  # The merge changes the header and group column names,
    # so they have to be changed back
    return out

def bin_overlap(data, left, right):
    # left_goneThrough = []  # X value of array that have already been . Ensures that there is only one of
    # each bin of bins (eg. [F][G1, G2, Gn])
    left_binsOfBins = []  # X value. Contains groups/families that have overlapping families/groups
    right_containedBins = []  # Y value. Contains families/groups that overlap with groups/families
    numSequences = []
    for i in data.index:
        leftFam = data[left + '_bins'].loc[i]
        rightGroup = data[right + '_bins'].loc[data[left + '_bins'] == leftFam].tolist()  # Finds
        # all instances of MEGARes groups that overlap with each CARD family
        if leftFam not in left_binsOfBins:
            left_binsOfBins.append(leftFam)
            right_containedBins.append(rightGroup)
            numSequences.append(len(rightGroup))
    df = {left: left_binsOfBins, right: right_containedBins, 'num_sequences': numSequences}
    return pd.DataFrame(df)


removeDuplicates = lambda x: list(dict.fromkeys(x))  # removes copies of the same string from a list


def num_bins(data, right):  # This needs to be done after duplicates have been removed, which is why it's not part of
    # bin_overlap
    numBins = []
    for i in data.index:
        numBins.append(len(data[right].loc[i]))
    df = pd.DataFrame({'num_bins': numBins})
    return pd.merge(data, df, left_index=True, right_index=True)

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

# Find percent of card sequences present in MEGARes and vice versa
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

# Get card-meg direction and remove duplicate bins
cardOverlap = bin_overlap(mergedDatabases, 'card', 'meg')
cardOverlap['meg'] = cardOverlap['meg'].apply(removeDuplicates)
# print(cardOverlap)

# Get meg-card direction and remove duplicate bins
megOverlap = bin_overlap(mergedDatabases, 'meg', 'card')
megOverlap['card'] = megOverlap['card'].apply(removeDuplicates)
# print(megOverlap)
# print(megOverlap.loc[megOverlap['card'].map(len) > 1])  # Finds MEG entries with more than one family

cardOverlap = num_bins(cardOverlap, 'meg')
megOverlap = num_bins(megOverlap, 'card')

print(cardOverlap)
print(megOverlap)

# print(megOverlap.loc[megOverlap['num_bins'] == 0])
# print(cardOverlap.loc[cardOverlap['num_bins']>1])

# DEBUG: Search for specific groups to test that bin_overlap is working properly
# searchgroup = 'AAC6-PRIME'  # MEG group to search for
# df2 = cardOverlap[[searchgroup in x for x in cardOverlap['meg']]]  # creates a dummy dataframe that holds all the
# # instances of that group that are present in cardOverlap
# print(cardOverlap.loc[df2.index])

# Print number of instances of bins of bins with i number of bins in them. Eg. 163 families have 1 group - 1: 163
def num_instances(data):
    numBins = []
    instanceList = []
    for i in range(0, 1000):
        numInstances = len(data[data['num_bins'] == i])
        if numInstances > 0:
            # print(i, end=' bin(s): ')
            # print(numInstances)
            numBins.append(i)
            instanceList.append(numInstances)
    return pd.DataFrame({'bins': numBins, 'instances': instanceList})

print("card instances: ")
print(num_instances(cardOverlap))
print("meg instances: ")
print(num_instances(megOverlap))

pd.DataFrame.to_csv(num_instances(cardOverlap),'card_num_bins.csv',index=False)