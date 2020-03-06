# - Do -PRIME notation conversion on CARD data
# - Connect each sequence to its CARD/MEGARes gene family - 2 dataframes
# Merge the dataframes along the sequence to produce a final dataframe:
# Dataframe: [Sequence, CARD gene fam, MEGARes group]
# - Figure out where the families and groups overlap - bar charts

# TODO: Find which groups are spread across multiple families and vice versa. Look at binsOfBins in card/megInstance
#  DataFrames.
# TODO: Remove Date stamping when translator is confirmed
# TODO: Get list of families in CARD that aren't in MEGARes and vice versa



import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from datetime import date

# Dataframe display options
# pd.options.display.width = 0
# pd.options.display.max_rows = 50
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
        g = lambda x: x[:x.index('_[')]  # removes species name. Can't just remove by space because some group names
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
        leftBin = data[left + '_bins'].loc[i]
        rightBin = data[right + '_bins'].loc[data[left + '_bins'] == leftBin].tolist()  # Finds
        # all instances of MEGARes groups that overlap with each CARD family
        if leftBin not in left_binsOfBins:
            left_binsOfBins.append(leftBin)
            right_containedBins.append(rightBin)
            numSequences.append(len(rightBin))
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


# Print number of instances of bins of bins with i number of bins in them. Eg. 163 families have 1 group - 1: 163
def num_instances(data, dbName):
    numBins = []
    instanceList = []
    binsOfBins = []
    for i in range(0, 1000):
        numInstances = len(data[data['num_bins'] == i])  # Go through the possible number of instances of each bin
        # within each bin of bins
        if numInstances > 0:
            # print(i, end=' bin(s): ')
            # print(numInstances)
            numBins.append(i)
            instanceList.append(numInstances)
            binsOfBins.append((data[dbName].loc[data['num_bins'] == i]).tolist())  # Adds a column containing all the
            # bins within that instance category (Eg. 10 bins, 1 instance, 16s rRNA methyltransferase (G1405)
    return pd.DataFrame({'num_bins': numBins, 'instances': instanceList, 'binsOfBins': binsOfBins})


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
numIn = 0
for i in cardData.index:
    if cardData['sequence'].loc[i] in megData['sequence'].tolist():
        numIn += 1

cardInMeg = (str(round(100 * numIn/len(cardData.index), 2)))
megInCard = (str(round(100 * numIn/len(megData.index), 2)))
print("Percent of CARD homolog model sequences in MEGARes: " + cardInMeg + "%")
print("Percent of MEGARes sequences in CARD homolog model: " + megInCard + "%")

#//region Find all group/family matches

mergedDatabases = pd.merge(left=cardData, right=megData, on='sequence')
mergedDatabases.rename(columns={'header_x': 'card_header', 'header_y': 'meg_header'}, inplace=True)

#//endregion

# Get card-meg direction and remove duplicate bins
cardOverlap = bin_overlap(mergedDatabases, 'card', 'meg')
cardOverlap['meg'] = cardOverlap['meg'].apply(removeDuplicates)

# Get meg-card direction and remove duplicate bins
megOverlap = bin_overlap(mergedDatabases, 'meg', 'card')
megOverlap['card'] = megOverlap['card'].apply(removeDuplicates)

# Get the number of bins contained in each bin of bins for each database
cardOverlap = num_bins(cardOverlap, 'meg')
megOverlap = num_bins(megOverlap, 'card')


# TODO: find bins in multiple bins of bins, then make this into a function and apply it to both CARD and MEG. Get
#  list of bins of bins that each bin falls into
#  - Create a separate DataFrame for it. Columns are bins of bins, rows are bins
#  - megInstances suggests that almost all groups fall into one family. Just analyze multi-overlap from MEG? No.
#  glycopeptide_resistance_gene_cluster;van_ligase has multiple families, as do APH3-PRIME and APH3-DPRIME

#//region Find bins that go into multiple bins of bins to figure out the differences in categorization between CARD
# and see which bins of bins are getting bins spread across them


# TODO: Translate this into a function so it can be applied in both directions

multiBinIndex = megOverlap.loc[megOverlap['num_bins'] > 1].index
megNames = megOverlap['meg'].loc[multiBinIndex]

collectCard = []
for i in megOverlap['card'].loc[megOverlap['num_bins'] > 1]:  # gets all the bins across all bins of bins. Only bins
    # in multiple bins of bins will show up more than once
    collectCard.extend(i)
dupCount = [x for x in collectCard if collectCard.count(x) > 1] # Only want to search for bins we KNOW are in
# multiple bins of bins, as bins that go into a single bin of bins will not be able to overlap
dupCount = removeDuplicates(dupCount)
# print(dupCount)

# df = pd.DataFrame(data=megOverlap['card'].loc[9], columns=[megOverlap['meg'].loc[9]])
# TODO: Rename df to be more descriptive
df = pd.DataFrame()
for i in multiBinIndex:  # Creates a dataframe that has the bins of bins as the column names and the bins as separate
    # rows. .loc[] has no use here. Have to search by using a for loop and megNames. This is done because Pandas has
    # no way to search an entire dataframe for a string, so you have to go through in a brute-force way
    columnData = megOverlap['card'].loc[i]
    # print(columnData)
    columnName = megOverlap['meg'].loc[i]
    df2 = pd.Series(data=columnData, name=columnName)
    # print(df2)
    df = pd.concat([df,df2],axis=1, ignore_index=True)
    # print(df)
df.columns = megNames
# print(df.loc['subclass_B3_LRA_beta-lactamase'])
print(df)

for duplicate in dupCount:
    dupCount = 0
    print("Bins of bins that contain " + duplicate + ": ", end=' ')
    binsOfBins = []
    for column in megNames:
        for row in df.index:
            if df[column].loc[row] == duplicate:
                dupCount += 1
                binsOfBins.append(df[column].name)
    print(binsOfBins)
    print("Number of bins of bins: " + str(dupCount))
    print('')
# print(df[megNames.tolist()[0]])

# df = pd.DataFrame(columns=megNames.tolist()) # Runs into index length mismatch issue
# for i in multiBinIndex:
#     df[megOverlap['meg'].loc[i]] = megOverlap['card'].loc[i]
#
#     print(df)

# print(megOverlap.loc[9])
# print(megOverlap['card'].loc[9][0:])
# print(cardOverlap.loc[cardOverlap['num_bins'] > 1])
# print(cardOverlap.loc[cardOverlap['card'] == 'OXA_beta-lactamase'])
exit()
# # print(len(megOverlap['card'].loc[megOverlap['num_bins'] > 1].tolist()))
# # print(len(megOverlap['meg'].loc[megOverlap['num_bins'] > 1].tolist()))
#
# dataIn = megOverlap['card'].loc[megOverlap['num_bins'] > 1].tolist()  # TODO: This doesn't work because it spreads
# # the data horizontally instead of vertically.
# print(dataIn)
# exit()
# test = pd.DataFrame(dataIn,
#                     columns=megOverlap['meg'].loc[megOverlap['num_bins'] > 1].tolist())
# print(test)
# print(megOverlap['card'].loc[megOverlap['meg'] == 'VANRA'])
# collectCard = []
# for i in megOverlap['card'].loc[megOverlap['num_bins'] > 1]:
#     collectCard.extend(i)
# # print(collectCard)
# for i in collectCard:  # DON'T USE .count IN A LOOP. PERFORMANCE ISSUE. USE Counter MODULE INSTEAD
#     if collectCard.count(i) > 0:
#         print(i)
#         print(megOverlap['card'].loc[megOverlap['card']])
#

# exit()
#//endregion


# # DEBUG: Search for specific groups to test that bin_overlap is working properly
# searchgroup = 'AAC6-PRIME'  # MEG group to search for
# df2 = cardOverlap[[searchgroup in x for x in cardOverlap['meg']]]  # creates a dummy dataframe that holds all the
# # instances of that group that are present in cardOverlap
# print(cardOverlap.loc[df2.index])

cardInstances = num_instances(cardOverlap, 'card')
megInstances = num_instances(megOverlap, 'meg')

print("card instances: ")
print(cardInstances)
print("meg instances: ")
print(megInstances)

print("EXIT EARLY - Just before to_csv")
exit()

filename = '_num_bins'

print("Writing CARD CSV...")
pd.DataFrame.to_csv(cardInstances,'card' + filename,index=False)
print("DONE")

print("Writing MEGARes CSV...")
pd.DataFrame.to_csv(megInstances,'meg' + filename,index=False)
print("DONE")