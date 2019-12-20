import csv

results = []
with open('aro_categories_index_reordered.tsv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    for row in reader: # each row is a list
        results.append(row)
# print(results[1])
# print("|".join(results[1]))
header = []

for i in range(0,len(results)):
    header.append("|".join(results[i]))

print(header)