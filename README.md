Given paper doi's, downloads structures from the Crystallography Open Database.
Will also download from the Cambridge Structure Database, if the user has a licensed installation of their Python API.
Individual molecules from each structure are outputted as .mol files.
A single molecule is operationally defined as a connected component of the bonding graph, using the bonding information in the .cif file from the database.

doi's are read from stdin.

Downloading molecules from one paper:

echo "[doi]" | python cod2mol.py

Downloading molecules from multiple papers, using a file containing each doi on a separate row:

python cod2mol.py < [name of file]

The doi will be extracted from each row, and rows without a doi will be skipped, so it's okay to have, for example, blank lines, or doi.org urls.

Information on the downloaded molecules is included in "output\_table.csv", which has columns "doi", "database", "entry", "molecule", and "filename".
If "output\_table.csv" already exists in the working directory, it will be appended to, and any doi already in it will not be redownloaded.
