import csv
import os
import os.path
import sys
import util


# Newline argument to prevent empty lines
# Following suggestion in this stackoverflow answer:
# https://stackoverflow.com/a/3348664/4434502
header_row = ["doi", "database", "entry", "filename"]
if os.path.isfile("output_table.csv") and \
    len([line for line in open("output_table.csv").readlines() if line != "\n"]) >= 2:
    prevtable = list(csv.reader(open("output_table.csv", "r", newline = "")))
    assert prevtable[0] == header_row
    doi_index = header_row.index("doi")
    doi_seen = set(row[doi_index] for row in prevtable[1:])
    entry_index = header_row.index("entry")
    entry_seen = set(row[entry_index] for row in prevtable[1:])
    output_table_file = open("output_table.csv", "a", newline = "")
    output_table = csv.writer(output_table_file)
else:
    doi_seen = set()
    entry_seen = set()
    output_table_file = open("output_table.csv", "w", newline = "")
    output_table = csv.writer(output_table_file)
    output_table.writerow(header_row)

for line in sys.stdin:
    doi = util.extract_doi(line)
    if doi in doi_seen or doi is None:
        continue
    doi_seen.add(doi)

    print("Looking up ids for doi {}".format(doi))
    print("Searching CSD for entries with doi {}".format(doi))
    entries = util.query_by_doi(doi)

    for entry in entries:
        util.save_entry_as_mol2(entry)
        row_dict = util.entry_to_row(entry)
        row_list = [row_dict[var] for var in header_row]
        output_table.writerow(row_list)

    if len(entries) == 0:
        row_dict = util.doi_to_empty_row(doi)
        row_list = [row_dict[var] for var in header_row]
        output_table.writerow(row_lit)

output_table_file.close()
