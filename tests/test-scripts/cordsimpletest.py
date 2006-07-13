# Richard Darst, 2005
#
# Test script for mp3.CordSimple
#

import numarray
import mp3
import mp3testmod

natype = numarray.Float32
frame0 = numarray.asarray( [[0, 0, 0],
                            [1, 1, 1],
                            [2, 2, 2],
                            [3, 3, 3]], natype )
frame1 = numarray.asarray( [[0, 0, 1],
                            [1, 1, 1],
                            [2, 2, 1],
                            [3, 3, 1]], natype )
frame2 = numarray.asarray( [[0, 0, 2],
                            [1, 1, 2],
                            [2, 2, 2],
                            [3, 3, 2]], natype )

C = mp3.CordSimple()
C.appendframe(frame0)
C.appendframe(frame1)
C.appendframe(frame2)
print "nframes should be 3:", C.nframes()
print "natoms should be 4:", C.natoms()
print "let's look at all three frames:"
# this uses the iterators
for frame in C.iterframes():
    print "frame ", C.framen(), ":"
    print frame

# test going back to the beginning:
C.zero_frame()

# Test writing it out.
C.writedcd("outputs/cordsimple_test1.dcd")
TF = mp3testmod.TestedFile("cordsimple_test1.dcd").test()
# Of course, you could also do frame by frame writing or whatever.


# test accessing our list 
print type(C.all_frames), len(C.all_frames)

