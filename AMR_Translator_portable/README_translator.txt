# The CARD-AMR++ Translator (CAT++)

CAT++ is a python script that converts CARD database and annotation files and produces AMR++-legible files. CAT++ takes CARD annotation information from the aro_index.tsv (downloaded from https://devcard.mcmaster.ca/tools/generatejson) to produce headers for sequences taken from nucleotide_fasta_protein_homolog_model.fasta, the protein homolog model of the CARD database.

# Usage

1. Place the Amalgamated_Translator_OneSource.py file into the same folder as the 'aro_index.tsv' and 'nucleotide_fasta_protein_homolog_model.fasta' files. The script must have the "Bio", "pandas", and "numpy" packages in the same folder in order to run
2. Run the script with the command 
	$ python Amalgamated_Translator_OneSource.py

Arguments:
-i: allow the user to specify an index file. Otherwise, the script will use any file with the name 'aro_index.tsv'
-d: allow the user to specify a database file. Otherwise, the script will use any file with the name 'nucleotide_fasta_protein_homolog_model.fasta'

# Output

NOTE: If you see a bunch of lines saying something like 
"Database entry on line 3099 Has been culled because its annotation was already assigned to another sequence. AMR++ only allows one header per sequence, and the conversion from gene to gene families results in a loss of granularity.
Database entry on line 3133 Has been culled because its annotation was already assigned to another sequence. AMR++ only allows one header per sequence, and the conversion from gene to gene families results in a loss of granularity."
This is normal. Nothing is broken. Read the step-by-step section below for an explanation.


The headers of each CARD database entry are converted from

	source|DNA Accession|Directionality(+/-)|Gene Start-Gene Stop|ARO:###|Gene [species of origin]

to

	DNA Accession|Type|Class|Mechanism|Group*

it places the annotation and database file in a folder named 'translations' within the same directory as the script.

Overlap_culled_DB.csv is the final file produced and includes the entries that were removed from the database because their DNA Accessions and gene families overlapped with those of other entries.

# Application to AMR++


If using AMR++ on linux, you MUST run dos2unix on the output files (translated database and annotation). It will not function otherwise.

If run on a singularity version of AMR++, the command to run AMR++ will resemble

	$ nextflow run main_AmrPlusPlus_v2.nf -profile singularity --reads "../*_{1,2}.fastq" --output CARD_mar18_db_sewage_sample_AMR++ --amr CARD_to_AMRplusplus_Database_2020_Mar_18.fasta --annotation CARD_to_AMRplusplus_Annotation_2020_Mar_18.csv

Some of the information to run AMR++ in this way is not currently (April 2020) available on the MEGARes Readme, so I wish you the best of luck if that command doesn't work.

# Step-by-Step and Error Messages

This section explains how the translator functions at the script level

1. The index file is read into pandas a dataframe and the database is read into a list

2. The columns necessary for the translation ['DNA Accession', 'Drug Class', 'Resistance Mechanism', 'AMR Gene Family', 'Model Name'] are read into the file. 'Model Name' represents the gene name.

3. Annotations that cannot be applied (because they are duplicates or, more critically, their gene and DNA Accessions overlap or they lack one of the necessary pieces of a header) are identified and removed from the annotation dataframe. They are collected (overlapRows) to have their corresponding database entries removed.

4. Database translation begins. DNA Accessions and genes are collected from the database to be compared to the annotation later. The header must be split apart to access the gene.

5. The genes are converted from ' notation to -PRIME notation to be in-line with MEGARes

6. The headers that indicate an error has occurred somewhere are defined

7. A pair of nested for loops sort headers (newHeaders) into the same order as the database and assign each slot in newHeaders the annotation that corresponds to the database entry in that same position.

8. granularityCulled collects all the headers that are present more than once and replaces the header with the granularityMessage flag, indicating that that header was removed. The loop goes through backwards so that only the first entry remains while all sequences that share its annotation are removed.

9. Any database entries found to overlap by DNA Accession and gene in their annotation, have a null value in their header, or lack an annotation are assigned an appropriate flag

10. The user is notified of which entries are removed and why. Line numbers are given for the DATABASE file ("Database entry on line ### ..."). copy the ARO # and ctrl-f for it in the index file to identify the problematic annotation. If these appear, it is ***NOT*** a sign that an error has occurred. The script will continue to function fine.

11. Any entries not assigned a header and not removed will trigger an error and notify the user.
	Error message: "ERROR: Some database entries are not being assigned headers, but are also not being culled. No files generated. They are:" [problematic entries]

A second check is performed to ensure that every entry has a unique header. Erroneous entries are displayed and the script will produce no files.
	Error message: "ERROR: Multiple copies of the same translated header are being created. No files have been generated. The following headers are being duplicated: " [problematic entries]

12. The final annotation file is written out with the current date appended to the end of its name

13. One final error check occurs to ensure that database entries are being properly removed when necessary. If any incorrect (overlapping, duplicate, null, etc) headers are still present in newHeaders, no database file is generated.
	"ERROR: culled headers are not being dropped from newHeaders and they are being sent to the database. "
          "Annotation file already generated, but no database file generated."

# Notes

* If the Protein Variant Model is added in the future, the |RequiresSNPConfirmation flag will need to be appended here. A section of commented-out code to do this can be found in the script (Amalgamated_Translator_OneSource.py) in the "Create AMR++-compliant header..." region. Note that that code will currently apply that flag to ALL entries put through the translator.

- legacy code dealing with "Protein Accession" is present from when these were needed to access gene family information. They remain in case future development of CARD causes more overlapping to the point where protein accessions are needed to sort.

# Future additions needed:

- "Database entry on line ### has been culled..." (granularityCulledMessage) needs to be cut down to a list of entries that have been cut. Don't spam the user.
- If Protein Accessions are needed again, don't forget to add a 'Protein Accession' column to newAnn.columns, newAnn = pd.concat( ... ), and aroCol to allow them to be imported and used.