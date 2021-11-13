import os
import sys
import re
import requests
import numpy as np
import csv
from slugify import slugify
from openbabel import pybel
from pymatgen.ext.cod import COD
from pymatgen.core import Structure
from pymatgen.core import Molecule
from pymatgen.io.cif import CifParser
from pymatgen.io.cif import str2float
from pymatgen.io.xyz import XYZ
from pymatgen.util.coord import pbc_shortest_vectors
import networkx as nx
from networkx.algorithms.traversal.breadth_first_search import bfs_edges
import pymysql

# Check if it is possible to search the Cambridge Structural Database, which
# requires that the user have a license
csd_available = True
try:
    import ccdc
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
metals = ["Ir", "Pt"]
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
    # readfile output:
    obabel_mol = list(pybel.readfile(format_a, outfile_temp))[0]
    obabel_mol.write(format_b, outfile_name)
    os.remove(outfile_temp)
    
output_table_file = open("output_table.csv", "w")
output_table = csv.writer(output_table_file)
output_table.writerow(["doi", "database", "entry", "molecule", "filename"])

for line in sys.stdin:
    doi = line.strip()
    # Retrieve Crystallography Open Database entry ID numbers for all entries
    # associated with the given doi by searching the database
    mysql_con = pymysql.connect(host=url, user="cod_reader", db="cod")
    mysql_cursor = mysql_con.cursor()
    
    print("Looking up ids for doi {}".format(doi))
    print("Searching COD for entries with doi {}".format(doi))
    structure_ids = query_executor(mysql_cursor, doi)

    for structure_id in structure_ids:
        print("Downloading structure {}".format(structure_id))
        # Download any cif files associated with the doi from the Crystallography Open
        # Database. 
        cif_text = \
            requests.get("http://{}/cod/{}.cif".format(url, structure_id)).text
        # Parse the cif file using pymatgen's cif parser
        parsed_cif = CifParser.from_string(cif_text)
        entries = list(parsed_cif._cif.data.keys())
        assert len(entries) == 1
        entry = entries[0]

        # Make a graph using the bonding information in the cif file
        bonding_graph = nx.empty_graph()
        try:
            for atom1, atom2 in zip(\
                parsed_cif._cif.data[entry].data["_geom_bond_atom_site_label_1"],
                parsed_cif._cif.data[entry].data["_geom_bond_atom_site_label_2"]):
                bonding_graph.add_edge(atom1, atom2)
        except KeyError:
            print("No bonding information in cif file")
            continue

        # Create a pymatgen structure object from the cif file
        data = parsed_cif._cif.data[entry].data
        atom_names = data["_atom_site_label"]
        lattice = parsed_cif.get_lattice(parsed_cif._cif.data[entry])
        species = data["_atom_site_type_symbol"]
        coords = np.array([[str2float(i) for i in data["_atom_site_fract_x"]],
            [str2float(i) for i in data["_atom_site_fract_y"]],
            [str2float(i) for i in data["_atom_site_fract_z"]]]).transpose()
        structure = Structure(lattice, species, coords,
            site_properties = {"label": atom_names})

        # Get the coordinates of the connected component around each metal atom,
        # in a cartesian coordinate system centered on the iridium
        metal_names = [i \
            for i in structure.site_properties["label"] \
            if re.match(metal_re, i) is not None]
        name2index = dict((label, i) \
            for i, label in enumerate(structure.site_properties["label"]))
        edge_lengths = dict()
        for metal_name in metal_names:
            print("Outputting molecule around metal atom {}".format(metal_name))
            component_coordinates = {metal_name: np.array([[[0., 0., 0.]]])}
            component_edges = bfs_edges(bonding_graph, metal_name)
            for edge in component_edges:
                displacement = pbc_shortest_vectors(structure.lattice,
                    structure.sites[name2index[edge[0]]].frac_coords,
                    structure.sites[name2index[edge[1]]].frac_coords)
                component_coordinates[edge[1]] = \
                    component_coordinates[edge[0]] + displacement
                edge_lengths[edge] = np.linalg.norm(displacement)
            element_symbols = [re.match("[A-Z][a-z]?", i).group(0) \
                for i in component_coordinates.keys()]
            molecule = Molecule(element_symbols,
                # The coordinate vectors have too many dimensions; flatten them
                list(i.flatten() for i in component_coordinates.values()))

            outfile_base = slugify("doi_{}_entry_{}_molecule_{}".\
                format(doi, structure_id, metal_name))
            XYZ(molecule).write_file(outfile_base + ".xyz")
            obabel_convert(outfile_base, "xyz", "mol")
            output_table.writerow([doi, "COD", structure_id, metal_name, outfile_base + ".mol"])

    # Try downloading from CSD if nothing was available from COD
    if len(structure_ids) == 0:
        print("Nothing found on COD.")
        if csd_available:
            print("Searching CSD for entries with doi {}".format(doi))
            csd_query = ccdc.search.TextNumericSearch()
            csd_query.add_doi(doi)
            with ccdc.io.MoleculeReader("CSD") as csd_reader:
                for hit in csd_query:
                    structure_id = hit.identifier
                    molecule_entry = csd_reader(structure_id)
                    for i, component in enumerate(molecule_entry.components):
                        outfile_base = slugify("doi_{}_entry_{}_molecule_{}".\
                            format(doi, structure_id, i))
                        with ccdc.io.MoleculeWriter(outfile_base + ".mol2") as writer:
                            writer.write(component)
                        obabel_convert(outfile_base, "mol2", "mol")
                        output_table.writerow([doi, "COD", structure_id,
                            i, outfile_base + ".mol"])

output_table_file.close()
