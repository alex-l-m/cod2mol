import csv
import os
import os.path
import sys
import re
import requests
import numpy as np
import pandas as pd
from slugify import slugify
from openbabel import pybel
from pymatgen.core import Structure
from pymatgen.core import Molecule
from pymatgen.io.cif import CifParser
from pymatgen.io.cif import str2float
from pymatgen.io.xyz import XYZ
from pymatgen.util.coord import pbc_shortest_vectors
import networkx as nx
from networkx.algorithms.traversal.breadth_first_search import bfs_edges
import pymysql
from rdkit import Chem

# Check if it is possible to search the Cambridge Structural Database, which
# requires that the user have a license
csd_available = True
try:
    import ccdc
    import ccdc.search
    import ccdc.io
except ModuleNotFoundError:
    print("ccdc module not installed, will only query COD")
    csd_available = False
except RuntimeError:
    # If the user has the ccdc module installed but does not have a license, it
    # will produce a RuntimeError
    print("ccdc module is installed but license may not be available, will only query COD")
    csd_available = False

# A regular expression for recognizing element symbols at the beginning of
# strings, which represent metals
metals = ["Ir"]
metal_re = "(" + "|".join(metals) + ")"

url = "www.crystallography.net"

def query_executor(cursor, doi):
    sql = 'select file from data where DOI like %s'
    cursor.execute(sql, (doi))

    structure_ids = [i[0] for i in list(cursor.fetchall())]
    return structure_ids

def obabel_convert(outfile_base, format_a, format_b):
    # I don't know how to output in mol, so save as xyz
    # then convert to mol
    outfile_temp = outfile_base + "." + format_a
    outfile_name = outfile_base + "." + format_b
    # Iterating over the "readfile" output seems to cause OSErrors on
    # Windows due to the file not being closed. Therefore, read the
    # whole iterator as a list without storing a reference to the
    # readfile output
    obabel_mols = list(pybel.readfile(format_a, outfile_temp))
    # 1169841 on CSD seems to have no coordinates in the entry. In this case,
    # the list will be empty
    try:
        obabel_mol = obabel_mols[0]
    except IndexError:
        print("OpenBabel did not find a molecule")
        return False
    obabel_mol.write(format_b, outfile_name, overwrite = True)
    os.remove(outfile_temp)
    return True

# Newline argument to prevent empty lines
# Following suggestion in this stackoverflow answer:
# https://stackoverflow.com/a/3348664/4434502
header_row = ["doi", "database", "entry", "filename"]
if os.path.isfile("output_table.csv"):
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
    output_table.writerow()

for line in sys.stdin:
    doi = util.extract_doi(line)
    if doi in doi_seen:
        continue
    doi_seen.add(doi)

    print("Looking up ids for doi {}".format(doi))
    print("Searching CSD for entries with doi {}".format(doi))
    entries = util.query_by_doi(doi)

    for entry in entries:
        save_entry_as_mol2(entry)
        row = entry_to_row(entry)
        output_table.writerow(row)

    if len(entries) == 0:
        row = doi_to_empty_row(doi)
        output_table.writerow(row)

output_table_file.close()
