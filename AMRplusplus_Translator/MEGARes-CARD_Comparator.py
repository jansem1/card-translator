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
pd.options.display.width = 0
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



def find_spread(data, binsOfBinsColumn, binsColumn, binsOfBinsType='bins of bins'):

    multiBinIndex = data.loc[data['num_bins'] > 1].index
    binsOfBinsNames = data[binsOfBinsColumn].loc[multiBinIndex]

    collectBins = []
    for bin in data[binsColumn].loc[data['num_bins'] > 1]:  # gets all the bins across all bins of bins. Only bins
        # spread across multiple bins of bins will show up more than once
        collectBins.extend(bin)
    dupCount = [bin for bin in collectBins if collectBins.count(bin) > 1]  # Only want to search for bins we KNOW are in
    # multiple bins of bins, as bins that go into a single bin of bins will not be able to overlap
    dupCount = removeDuplicates(dupCount) # remove duplicates so we only go once through each bin of bins across which bins
    # are spread
    # print(dupCount)
    binsInDataFrame = pd.DataFrame()
    for row in multiBinIndex:  # Creates a dataframe that has the bins of bins as the column names and the bins as
        # separate rows. Have to search by using a for loop and binsOfBinsNames. This is
        # done because Pandas has no way to search an entire dataframe for a string, so you have to go through in a
        # brute-force way.
        columnData = data[binsColumn].loc[row]
        df2 = pd.Series(data=columnData)
        binsInDataFrame = pd.concat([binsInDataFrame,df2],axis=1, ignore_index=True)
    binsInDataFrame.columns = binsOfBinsNames
    # print(binsInDataFrame.loc['subclass_B3_LRA_beta-lactamase'])
    # print(binsInDataFrame)
    x = []
    y = []
    for duplicate in dupCount:
        dupNumber = 0
        binsOfBin = []  # Bins of bins across which the current bin (duplicate) is spread across
        for column in binsOfBinsNames:
            for row in binsInDataFrame.index:
                if binsInDataFrame[column].loc[row] == duplicate:
                    dupNumber += 1
                    binsOfBin.append(binsInDataFrame[column].name)
        print(duplicate + " (" + binsColumn + ")" + " is spread across these " + binsOfBinsType + ": ", end=' ')
        print(binsOfBin)
        print("Number of " + binsOfBinsType + ": " + str(dupNumber) + '\n')
        x.append(duplicate)
        y.append(binsOfBin)
    return pd.DataFrame({binsColumn:x, 'spread_' + binsOfBinsType: y})


#//endregion

#//region Comparison

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

#//region Find bins that go into multiple bins of bins to figure out the differences in categorization between CARD
# and see which bins of bins are getting bins spread across them

megSpread = find_spread(megOverlap,'meg','card','groups')
cardSpread = find_spread(cardOverlap,'card','meg','families')

print(megSpread)
print(cardSpread)

#
# print("EXIT EARLY - Just before to_csv")
# exit()

instanceFilename = '_num_bins.csv'
spreadFilename = '_spread.csv'

print("Writing CARD CSVs...")
pd.DataFrame.to_csv(cardInstances,'card' + instanceFilename,index=False)
print("Instances DONE")
pd.DataFrame.to_csv(cardSpread, 'card' + spreadFilename, index=False)
print("Spread DONE")

print("Writing MEGARes CSV...")
pd.DataFrame.to_csv(megInstances,'meg' + instanceFilename,index=False)
print("Instances DONE")
pd.DataFrame.to_csv(megSpread, 'meg' + spreadFilename, index=False)
print("Spread DONE")
#//endregion
#//endregion