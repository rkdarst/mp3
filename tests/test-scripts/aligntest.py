import mp3

inputfile = "inputs/water1000.dcd"
psffile = "inputs/water1000.psf"
#inputfile = "/mnt/pete29/home/richard/lonestar/apomb_310/masterDCD/apomb_eq-32-36_WRAP.dcd"
#inputfile = "/home/richard/tmpdcd.dcd"
#psffile = "apomb_SPCE100ion.psf"

# Load up the system (labels data)
S = mp3.System()
print "reading psf"
S.labels.getfrompsf( psffile )
print "done reading psf"

# First, we've got to find which atoms we want to align with.
print "finding atomlist"
atomlist = list(S.labels.findatoms(atomtype="OSPC"))
#print atomlist
print "done finding atomlist"



# First off, a basic test.  Align by some chosen atoms.  Have a hook
# called after each frame.
if True:
    # First, create a cord to read from the DCD
    C = mp3.CordDCD(dcd=inputfile )
    #C = mp3.CordDCD(dcd="/mnt/pete29/home/clopez/research/apombC/310Kapomb2C/voronoi/frame100-200_slice.dcd" )
    C = mp3.CordDCD(dcd="inputs/water1000.dcd" )
    atomlist = range(2464)
    # This creates the alignment object.
    C = mp3.CordAlign(cord=C,
                      atomlist=atomlist)
    C.writedcd("tmp-aligntest.dcd")
    C.zero_frame()
    # This is jus a way to let us do stuff after each alignment.  It's
    # not very often used.
    def nextframe_end_hook(self):
        print self.error,
        print self.iterations
    C.nextframe_end_hook = nextframe_end_hook

    # Do the reading.
    print "reading first frame"
    C.nextframe()
    for i in range(8):
        C.nextframe()
    saved_frame = C.nextframe()
    print "Test 1 completed"

    

# Test 2.  Align by all atoms and take the average frame
if True:
    C = mp3.CordDCD(dcd=inputfile)
    C = mp3.CordAlign(cord=C,
                      saveaverage= True)
    C.nextframe()
    C.averageframe()
    C.read_n_frames(8)
    C.averageframe()
    print "Test 2 completed"



# Test 3.  Align by all atoms, aligining to a previously-known frame.
if True:
    C = mp3.CordDCD(dcd=inputfile )
    C = mp3.CordAlign(cord=C,
                      aligntoframe=saved_frame)
    C.nextframe()
    C.read_n_frames(8)
    C.nextframe()
    print "Test 3 completed"



# Test 4.  Align by all atoms, aligining to a previously-known frame.
if True:
    C = mp3.CordDCD(dcd=inputfile )
    C = mp3.CordAlign(cord=C,
                      atomlist=atomlist,
                      aligntoframe=saved_frame[atomlist])
    C.nextframe()
    C.read_n_frames(8)
    C.nextframe()
    print "Test 4 completed"
