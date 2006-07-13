#!/usr/bin/env python

import mp3

# Open our input DCD.

C = mp3.CordDCD(dcd="inputs/water1000.dcd")
print "created DCD object, natoms:", C.natoms()


# So, begin our step-by-step writing.  The first argument is the
#  file we are writing to.  The optional keyword argument "nframes"
#  is because the original DCD has 10 frames, but we are only
#  writing 5 of them.  The number of frames is encoded in the header,
#  so we have to set it properly before we start.
#    NOTE:  a more advanced method is planned.

C.writedcd_start("outputs/water1000_fbf.dcd", nframes=5)

for i in range(5):
    C.nextframe()         # You have to call nextframe() each time.

    # Do something with the frame:
    print "first five atoms of frame %s: "%i*2
    print C.frame()[0:5]

    
    C.writedcd_nextframe()   # Write the current frame.


C.writedcd_stop()            # We have to close the file
                             #  (at least you probably want to)
