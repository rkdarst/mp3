
import mp3

C = mp3.CordPDB("inputs/apombC_SPCE100ion.pdb")
print C.framen()
print C.natoms()

C.nextframe()

print C.frame().shape

C.natoms()
C.nframes()
C.frame()
C.firsttstep()
C.dcdfreq()
C.tstep_size()
C.block_a()
C.block_b()
C.charm_v()
C.title()


# test writing a PDB
S = mp3.System("inputs/water1000.psf",
               cord=mp3.CordDCD("inputs/water1000.dcd"))
S.writepdbseries("outputs/pdbs/water1000--")

