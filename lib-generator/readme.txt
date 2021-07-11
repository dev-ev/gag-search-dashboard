The tab delimited files are a human-readable starting point for library generation,
tab files contain the following columns:
1) unique name of the glycan
2) number of sulfates
3) number of hexuronic acid residues
4) number of NeuAc residues
5) neutral monoisotopic mass of the molecule
6) minimal retention time in minutes
7) maximal retention time in minutes

The script reads the tabulated input and produces the libraries of the theorethical GAG forms,
including different charge states and dibutylamine adducts. The databases are saved as JSON files.

The examples of the JSON libraries for chondroitin/dermatan sulfates (CS/DS) and heparan sulfate (HS)
can be found in the zip file in this folder.
