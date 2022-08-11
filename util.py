import re
import ccdc
import ccdc.search
import ccdc.io

# DOI regex from:
# https://www.crossref.org/blog/dois-and-matching-regular-expressions/
# The first regex given fails on this one:
# 10.1002/1521-4095(20020116)14:2<147::AID-ADMA147>3.0.CO;2-3
# Therefore, using the second one, and then trying the first (modern) one,
# modified to be case insensitive
old_doi_regex = r"10.1002/[^\s]+"
new_doi_regex = r"10.\d{4,9}/[-._;()/:A-Za-z0-9]+"

# Check that an input is a string containing a doi
def validate_doi(string):
    if not isinstance(string, str):
        raise ValueError("DOI is not a string")
    if re.match(old_doi_regex + "$", string) is None \
            and re.match(new_doi_regex + "$", string) is None:
        raise ValueError(f"DOI {string} is not valid")

def extract_doi(string):
    old_doi_match = re.search(old_doi_regex, string)
    if old_doi_match is not None:
        return old_doi_match.group(0)
    new_doi_match = re.search(new_doi_regex, string)
    if new_doi_match is not None:
        return new_doi_match.group(0)
    return None

ccdc_id_regex = "[A-Za-z][A-Za-z][A-Za-z][A-Za-z][A-Za-z][A-Za-z]"

def validate_ccdc_id(string):
    if not isinstance(string, str):
        raise ValueError("DOI is not a string")
    if re.match(ccdc_id_regex + "$", string) is None:
        raise ValueError(f"CCDC ID {string} is not valid")

def extract_ccdc_id(string):
    ccdc_id_match = re.search(ccdc_id_regex, string)
    if ccdc_id_match is not None:
        return ccdc_id_match.group(0)
    return None

def query_by_doi(doi):
    validate_doi(doi)
    csd_query = ccdc.search.TextNumericSearch()
    csd_query.add_doi(doi)
    # List of "SearchHit" objects
    csd_results = list(csd_query.search())
    entries = [hit.entry for hit in csd_results]
    return entries

def query_by_id(ccdc_id):
    validate_ccdc_id(ccdc_id)

    with ccdc.io.EntryReader("CSD") as csd_entry_reader:
        entry = csd_entry_reader.entry(ccdc_id)

    return entry

def save_entry_as_cif(entry):
    # Doesn't seem to be the same as the cif file you get from the website
    # Probably not going to use this
    filename = entry.identifier + ".cif"
    cif_file_text = entry.to_string("cif")
    open(filename, "w").write(cif_file_text)

def save_mol_as_mol2(entry, mol):
    filename = entry.identifier + ".mol2"
    with ccdc.io.MoleculeWriter(filename) as writer:
        writer.write(mol)

def entry_to_row(entry, mol):
    row = dict({\
        "doi": entry.publication.doi,
        "database": "CSD",
        "entry": entry.identifier,
        "filename": entry.identifier + ".mol2",
        "crystal_formula": entry.formula,
        "molecule_formula": mol.formula,
        "formal_charge": mol.formal_charge})
    return row

def doi_to_empty_row(doi):
    row = dict({\
        "doi": doi,
        "database": None,
        "entry": None,
        "filename": None,
        "crystal_formula": None,
        "molecule_formula": None,
        "formal_charge": None})
    return row

def entry_to_empty_row(entry):
    row = dict({\
        "doi": None,
        "database": "CSD",
        "entry": entry,
        "filename": None,
        "crystal_formula": None,
        "molecule_formula": None,
        "formal_charge": None})
    return row

def entry_to_components(entry):
    return entry.molecule.components

def molecule_element_count(mol, symbol):
    return len([atom for atom in mol.atoms if atom.atomic_symbol == symbol])
