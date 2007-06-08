#!/usr/bin/env python


import mp3
import mp3testmod

# Create DCD object
C = mp3.CordDCD(dcd="inputs/water1000.dcd")
assert C._fo == None

# Some important properties (look at mp3testmod.py to see)
mp3testmod.cord_attributes_test(C)

# write DCD
C.writedcd("outputs/water1000.dcd")

# DONE WITH TEST

# verify it
mp3testmod.TestedFile("water1000.dcd")


C.zero_frame()
assert C._fo == None
for i in range(C.nframes()):
    C.nextframe()
assert C._fo == None

C.zero_frame()
C.read_n_frames(C.nframes())
assert C._fo == None
