#!/usr/bin/env python

#
# Labels test script.  Most of this comes from the lesmod test script.
#

import sys
import sets
import shelve

import mp3


print "loading pdb... ", 
S = mp3.System(pdb="inputs/les/apombC_SPCE100ion.pdb")
print "pdb loaded"

# the replicated residues are: 1-5, 19-107, 110, 120-126, 138, 140-153
# let's find which atoms are replicated. (see the python help for sets)
print "searching for first atomset... ",
atoms = sets.Set()
atoms |= S.findatoms(resnum=(1,5))
print "first atomset found"
atoms |= S.findatoms(resnum=(19,107))
atoms |= S.findatoms(resnum=110)
atoms |= S.findatoms(resnum=(120,126))
atoms |= S.findatoms(resnum=138)
atoms |= S.findatoms(resnum=(140,153))
atoms &= S.findatoms(segname="APO")
atomlist = list(atoms)
atomlist.sort()
print atomlist


# we have to be sure that findrange is still backwards-compatable

S.findrange(resnum=(19,107))


   
print "atomset created"
