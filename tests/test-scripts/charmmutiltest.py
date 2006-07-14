# Richard Darst, 2006

import mp3.log
import mp3.charmmutil

mp3.log.set_level("debug")
filename = "inputs/charmm/top_all27_prot_lipid.inp"
Masses = mp3.charmmutil.getMasses(filename)
print Masses["HE1"]
print Masses["HE1"]["symbol"]
print Masses["HE1"]["number"]
print Masses["HE1"]["mass"]

Masses = mp3.charmmutil.getMasses(filename, sortby="number")
Masses = mp3.charmmutil.getMasses(filename, sortby="symbol")
Masses = mp3.charmmutil.getMasses(filename, sortby="mass")


Residues = mp3.charmmutil.getResidues(filename)
print Residues["ALA"]
print Residues["ALA"]["charge"]
print Residues["ALA"]["atomCharge"]["CA"]
print Residues["ALA"]["atomNamesByType"]["CT1"]
print Residues["ALA"]["atomNames"]
print Residues["ALA"]["atomTypes"]
print Residues["ALA"]["groups"]

#from rkddp.interact import interact; interact()

