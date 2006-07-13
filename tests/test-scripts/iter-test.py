
import mp3

# Create DCD object
C = mp3.CordDCD(dcd="inputs/water1000.dcd")
#print C.nframes()

print "testing itercords"
for C in C.itercords():
    print C.framen()
print


C = mp3.CordDCD(dcd="inputs/water1000.dcd")
print "testing iterframes"
for frame in C.iterframes():
    # since C refers to the same object as we go through, we can
    # continue to use C throught the iterator, even though the
    # iterator is over frames.
    print C.framen()
    print frame[0] 
print
