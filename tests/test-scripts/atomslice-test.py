#!/usr/bin/env python


import mp3
import mp3testmod

# Create DCD object
C = mp3.CordDCD(dcd="inputs/water1000.dcd")

# Some important properties (look at mp3testmod.py to see)
C= mp3.CordAtomSlice(cord=C, atomlist=range(500))

mp3testmod.cord_attributes_test(C)
print C.nframes()

for i in range(C.nframes()):
    C.nextframe()

