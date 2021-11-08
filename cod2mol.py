import sys
import re
import requests
import numpy as np
from slugify import slugify
from pymatgen.ext.cod import COD
from pymatgen.core import Structure
from pymatgen.core import Molecule
from pymatgen.io.cif import CifParser
from pymatgen.io.cif import str2float
from pymatgen.io.xyz import XYZ
from pymatgen.util.coord import pbc_shortest_vectors
import networkx as nx
from networkx.algorithms.traversal.breadth_first_search import bfs_edges

# A regular expression for recognizing element symbols at the beginning of
# strings, which represent metals 
metals = ["Ir", "Pt"]
metal_re = "(" + "|".join(metals) + ")"

for line in sys.stdin:
    doi = line.strip()
    # Retrieve Crystallography Open Database entry ID numbers for all entries
    # associated with the given doi by searching the database
    downloader = COD()
    print("Looking up ids for doi {}".format(doi))
    sql_request = "select file from data where DOI like '{}'".format(doi)
    structure_ids = [int(i) for i in downloader.query(sql_request).split()[1:]]

    for structure_id in structure_ids:
        print("Downloading structure {}".format(structure_id))
        # Download any cif files associated with the doi from the Crystallography Open
        # Database. 
        cif_text = \
            requests.get("http://{}/cod/{}.cif".format(downloader.url, structure_id)).text
        # Parse the cif file using pymatgen's cif parser
        parsed_cif = CifParser.from_string(cif_text)
        entries = list(parsed_cif._cif.data.keys())
        assert len(entries) == 1
        entry = entries[0]

        # Make a graph using the bonding information in the cif file
        bonding_graph = nx.empty_graph()
        for atom1, atom2 in zip(\
            parsed_cif._cif.data[entry].data["_geom_bond_atom_site_label_1"],
            parsed_cif._cif.data[entry].data["_geom_bond_atom_site_label_2"]):
            bonding_graph.add_edge(atom1, atom2)

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
            element_symbols = [re.match("[A-Z][a-z]?", i).group(0) for i in component_coordinates.keys()]
            molecule = Molecule(element_symbols,
                # The coordinate vectors have too many dimensions; flatten them
                list(i.flatten() for i in component_coordinates.values()))

            outfile_name = slugify("doi_{}_entry_{}_molecule_{}".\
                format(doi, structure_id, metal_name)) + ".xyz"

            # Write a single molecule from the crystal as an output xyz file
            XYZ(molecule).write_file(outfile_name)
