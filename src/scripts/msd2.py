

import hotshot
import shelve
import numarray
import mp3
import Gnuplot

def main():

    # Put what you want here.
    dcd_to_msd = "/home/richard/research/mp3/tests/water1000.dcd"
    psf_name = "/home/richard/research/mp3/tests/water1000.psf"
    
    # Set up the systems, so that you can find the atoms you want below.

    cord = mp3.dcd()
    cord.setinput(dcd_to_msd)
    s = mp3.system(cord=cord, psf=psf_name)

    # Find which atoms you want right here.
    # First, we find all atoms in the protein, and all in the residues we
    # want.
    
    residue_low = 10
    residue_high = 20
    residues = [i for i in xrange(s.labels.natoms) if
                (s.labels.data.field("resnum") >= residue_low and
                 s.labels.data.field("resnum") <= residue_high)    ]
    protein = s.labels.findatoms(name="SPCE", field="resname")
    
    # Find alpha carbons, add them all to a set.  This is what we align by
    
    alpha_carbons = s.labels.findatoms(name="CA", field="atomtype")
    atomset = sets.Set()
    atomset |= protein  # union update
    atomset &= alpha_carbons  # intersection update
    atomset &= residues
    atoms_to_align_by = list(atomset)

    # Find all atoms in the residues.  Then find all hydrogens, and remove 
    # the hydrogens from the set
    
    # residues has been defined above
    hydrogens = s.labels.findatoms(name="H", field="atomtype")
    atomset = sets.Set()
    atomset |= protein
    atomset &= residues
    atomset -= hydrogens
    atoms_to_msd = list(atomset)
    
    # Set up our aligned cord, now

    #atomlist = range(cord.natoms)
    #atoms_to_align_by = atomlist
    #atoms_to_msd = atomlist
    
    aligned_cord = mp3.aligncord()
    aligned_cord.setcord(s.cord)
    aligned_cord.setatoms( atoms_to_align_by )
    aligned_cord.init()
    s.cord = aligned_cord

    msd = Msder(cord=s.cord, window=5, atomlist=atoms_to_msd )
    

    print "beginning correlation"
    msd.do_n_msd(5)
    

    print msd.msd()
    

class Msder:

    def __init__(self, cord, window, atomlist):

        # save some class variables
        self.cord = cord
        self.window = window
        self.atomlist = atomlist

        #make the storage we need for our stuff
        self.framelist = []
        self._msd_sum = [0] * self.window      # this stores the sum of all msd's...
        self._count_sum = [0] * self.window    # this stores the number of msd's at each distance

        # Load up what's in the window
        for i in range(self.window):
            self._append_frame()

## stuff for you to use

    def finish_msd_window(self):
        for i in xrange(self, len(self.framelist)):
            del self.framelist[0]
            self._one_msdrun()

    def do_n_msd(self,n=1):
        for i in xrange(n):
            self._advance_frame()
            self._one_msdrun()

    def do_all_msd(self):
        """Does as many MSD as this instance will allow.

        """
        n = self.cord.nframes - (self.cord.framen+1)
        do_all_msd(n)

    def msd(self):
        return numarray.divide(self._msd_sum, self._count_sum)

    def write_msdlist(self, filename="msd.dat"):
        fo = file(filename,"w")
        for num in self.msd():
            fo.write("%s\n"%num)


## private stuff


    def _append_frame(self):
        print None
        self.framelist.append(self.cord.nextframe()[self.atomlist].copy())

    def _advance_frame(self):
        """Advances the contents of self.framelist

        pulls one 
        """
        self._append_frame()
        del(self.framelist[0])

    def _one_msdrun(self):
        """Runs msd for self.framelist, comparing the first frame with each other frame
        """
        self._msd_sum[0] += 1
        for i in range(1, len(self.framelist)):     # go through our saved frames
            self._msd_sum[i] += mp3.rmsd(self.framelist[0], self.framelist[i])**2  # find rmsd
            self._count_sum[i] += 1

        

    # def shelve_msd(self, filename="state.pickle"):
    #     self.calc_msd()
    #     shelf = shelve.open(filename)
    #     shelf["msdlist"] = self.msdlist
#   #      shelf["msdobj"] = self
    #     shelf.close()

    


main()
