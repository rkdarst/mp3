import sys

import mp3
import mp3.parsecm3d
import mp3.cordcm3d

import mp3.log
mp3.log.set_level("debug")

#from hotshot import Profile
#profile = Profile("profile")

#
# Test reading set files.
#

#profile.start()
S = mp3.System()
S.readcm3dset("inputs/cm3dfiles/apo_wat.set")
#profile.stop()
print S.labels.data[:10]
# sys.exit()

# Set the protein atom names.  We read them from another psf file we
# have, and we know that there are 2464 protein atoms.
Spsf = mp3.System(psf="inputs/cm3dfiles/apomb_works.psf")
print Spsf.data["atomname"]
S.data['atomname'][0:2464] = Spsf.data['atomname']
# Set the water atomnames and resnames.  We assume that all remaining
# atoms are water.
nWat = len(S.data['atomname'][2464:])/3
S.data['atomname'][2464:] = ['OH','1H','2H'] * nWat
S.data['resname'][2464:] = ['HOH'] * nWat * 3

#
# Test trajectory files.
#

S.setcord(mp3.cordcm3d.CordCM3D("inputs/cm3dfiles/apo-wat-10.confp"))

S.cord.nextframe()

# Test writing a PBD.  This uses the current coordinates and all of
# the labels we have loaded.  Any other PBD writing functions should
# work as well.
S.writepdb("test.pdb")
#from rkddp.interact import interact ; interact()
#sys.exit()


S.cord.read_n_frames(8)
S.cord.nextframe()
S.cord.zero_frame()
assert S.cord.framen() == -1

print S.cord.nframes()
S.cord.calcnframes()
print S.cord.nframes()
assert S.cord.nframes() == 10

S.cord.nextframe()
S.cord.nextframe()
S.cord.nextframe()
S.cord.nextframe()
S.cord.nextframe()
S.cord.nextframe()
S.cord.nextframe()
S.cord.nextframe()

S.cord.calcnframes()
print "the next line should give an error and only write two frames"
S.cord.write_cm3d("outputs/cm3d-apo-wat.confp", maxframes=11)



# Test the saveFD and __init__(ensure_nframes=True) functionality.
print "\nbeginning saveFD and __init__(ensure_nframes=True) tests:"
C = mp3.cordcm3d.CordCM3D("inputs/cm3dfiles/apo-wat-10.confp",
                          ensure_nframes=True, saveFD=True)
print C.nframes()
for x in range(10):
    C.nextframe()
print "done"



#from rkddp.interact import interact ; interact()
