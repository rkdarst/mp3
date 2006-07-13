#!/usr/bin/env python

import sys
import sets
import mp3
import shelve
shelf = shelve.open("inputs/les/lesatomlist.shelve")

use_stored = True
if use_stored == True:
    atomlist = shelf["atomlist"]
else:
    S = mp3.System(pdb="inputs/les/apombC_SPCE100ion.pdb")
    print "pdb loaded"
    
    # the replicated residues are: 1-5, 19-107, 110, 120-126, 138, 140-153
    # let's find which atoms are replicated. (see the python help for sets)
    
    atoms = sets.Set()
    atoms |= S.findrange(resnum=(1,5))
    print "first atomset found"
    atoms |= S.findrange(resnum=(19,107))
    print atoms
    atoms |= S.findatoms(resnum=110)
    atoms |= S.findrange(resnum=(120,126))
    atoms |= S.findatoms(resnum=138)
    atoms |= S.findrange(resnum=(140,153))
    atoms &= S.findatoms(segname="APO")
    atomlist = list(atoms)
    atomlist.sort()
    print atomlist
    #shelf["atomlist"] = atomlist
   
print "atomset created"


# Read in the LES DCD.  This should be like normal.

C = mp3.CordDCD(dcd="inputs/les/apombC_310KLESnc_pre1.dcd")
print "created DCD object, natoms:", C.natoms()


# Now, we want to Un-LES it.  mp3.CordLESMod is a wrapper that
#  un-LES'es it.  
#    - Read coordinates from C
#    - numrep=5 specifies each replicated atom is replicated 5 times
#    - atomlist is the
#    - whichframe=-1  tells it that you want the average (see help (mp3.CordLESMod))

C = mp3.CordLESMod(cord=C,                     
                   numrep=5,
                   #atomlist=atomlist,
                   leslist = shelf["leslist"],
                   whichframe = -1)
print "created LES object, natoms:", C.natoms()


#shelf["leslist"] = C._leslist
#shelf.close()

print C._leslist

print sum(C._leslist)   # both of these should give the same thing.
print C.natoms()        #

print "nframes:", C.nframes()


#C.writedcd("tep.dcd")
C.nextframe()

C.rmsd()
C.variance()
C.mean()
C.frames()
C.frame()
print C.rep_rmsd(1, atomlist)



#print "nextframe completed successuflly"
#for i in range(100):
#    C.nextframe()
#    print i
    
