import numarray
import mp3

s = mp3.system()
s.cord = mp3.dcd()
s.cord.setinput("water1000.dcd")
s.labels.getfrompsf("water1000.psf")

s.cord.nextframe()

# return a list of protein atoms:
 # make a list containing atom numbers and properties

atoms = [ i for i in xrange(0,s.cord.natoms) if s.labels.data.field('atomtype')[i] == "OSPC" ]

goodatoms = s.cord.frame[atoms]
#print goodatoms

minatoms = goodatoms.argmin(axis=0)
maxatoms = goodatoms.argmax(axis=0)

boxmin = numarray.asarray( [ goodatoms[minatoms[i],i] for i in (0,1,2) ] )
boxmax = numarray.asarray( [ goodatoms[maxatoms[i],i] for i in (0,1,2) ] )

numarray.add(boxmin, -10, boxmin )
numarray.add(boxmax,  10, boxmax )

print boxmin, boxmax
