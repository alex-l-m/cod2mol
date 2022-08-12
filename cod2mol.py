import csv
import os
import os.path
import sys
import util


# Newline argument to prevent empty lines
# Following suggestion in this stackoverflow answer:
# https://stackoverflow.com/a/3348664/4434502
header_row = ["doi", "database", "entry", "deposition_number", "filename", "crystal_formula", "molecule_formula", "formal_charge"]
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

    print("Searching CSD for entries with doi {}".format(doi))
    entries = util.query_by_doi(doi)

    for entry in entries:
        components = util.entry_to_components(entry)
        components_with_one_ir = [component for component in components \
                if util.molecule_element_count(component, "Ir") == 1]
        if len(components_with_one_ir) == 0:
            print("No components with one iridium")
            row_dict = util.doi_to_empty_row(doi)
            row_list = [row_dict[var] for var in header_row]
            output_table.writerow(row_list)
        else:
            print("Found a component with one iridium")
            component_to_write = components_with_one_ir[0]
            util.save_mol_as_mol2(entry, component_to_write)
            row_dict = util.entry_to_row(entry, component_to_write)
            row_list = [row_dict[var] for var in header_row]
            output_table.writerow(row_list)

    if len(entries) == 0:
        print("No results from search")
        row_dict = util.doi_to_empty_row(doi)
        row_list = [row_dict[var] for var in header_row]
        output_table.writerow(row_list)

output_table_file.close()
