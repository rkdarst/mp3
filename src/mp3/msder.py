

import numarray
import mp3.functions

class Msder:
    """Perform mean-square calculations on coordinate data.
    """
    def __init__(self, cord, window, atomlist=None):
        """Set coordinate set and window ; load first `window' frames.
        """
        # if atomlist is None, then make it all atoms
        if atomlist == None:
            atomlist = range(cord.natoms())

        # save some class variables
        self.cord = cord
        self.window = window
        self.atomlist = atomlist
        self._degrees_of_freedom = 3           # used for diffusion constant

        #make the storage we need for our stuff
        self.framelist = []
        self._msd_sum = [0] * self.window      # this stores the sum of all msd's...
        self._count_sum = [0] * self.window    # this stores the number of msd's at each distance

        # Load up what's in the window
        for i in range(self.window):
            self._append_frame()

## stuff for you to use

    def do_n_msd(self,n=1):
        """Advance the coordinates and do n msd's.
        """
        for i in xrange(n):
            self._advance_frame()
            self._one_msdrun()

    def do_all_msd(self):
        """Does as many MSD as this instance will allow.

        This does not "complete" the msd.  To get all possible data out
        (resulting in more samples at a shorter time distances), use
        finish_msd_window after this.

        """
        n = self.cord.nframes() - (self.cord.framen()+1)
        self.do_n_msd(n)

    def finish_msd_window(self):
        """`Complete' the msd of the system.
        """
        for i in xrange(len(self.framelist)):
            del self.framelist[0]
            self._one_msdrun()

    def msd(self):
        """Return a numarray containing the msd.
        """
        return numarray.divide(self._msd_sum, self._count_sum)

    def write_msdlist(self, filename="msd.dat"):
        """Write the msd list (separated by newlines) to a file.
        """
        fo = file(filename,"w")
        for num in self.msd():
            fo.write("%s\n"%num)


## private stuff


    def _append_frame(self):
        """Add a frame to the list of recent frames.
        """
        self.framelist.append(self.cord.nextframe()[self.atomlist].copy())
        #print "."

    def _advance_frame(self):
        """Advances the contents of self.framelist

        Appends a frame, then deletes the oldest one.
        """
        self._append_frame()
        del(self.framelist[0])

    def _one_msdrun(self):
        """Runs msd for self.framelist.

        Compare the first frame with each other frame, storing the results.
        """
        self._count_sum[0] += 1
        for i in range(1, len(self.framelist)):     # go through our saved frames
            self._msd_sum[i] += mp3.functions.rmsd(self.framelist[0], self.framelist[i])**2  # find rmsd
            self._count_sum[i] += 1


## to find the msd

    def set_distanceunit(self, unit):
        """Set what the distance unit on this system is.

        This should be the scale of the (cartesian) coordinate system
        being modeled.  A typical value would be "angstrom".  Hovewer,
        other forms would be acceptable, for example, "2 angstrom",
        if the distance of unity was in actuality 2 angstroms (note
        that this is a corner case designed to demonstrate the flexability)
        """
        self._distance_unit = unit

    def set_timeunit(self, unit):
        """Set what the time unit on this system is.

        This should be the interval between frames in this system.
        A typical value would be ".1 picosecond".
        """
        self._time_unit = unit
        
    def find_diffusion_constant(self, begin, end):
        """Find the diffusion constant.

        This method returns a tuple containing the slope and units
        of the slope.
        """
        msdlist = self.msd()
        timelist = range(len(msdlist))

        print "WARNING-- using simple two point slope approximation instead of least squares"
        slope = (msdlist[end]-msdlist[begin]) / (timelist[end]-timelist[begin])

        print "WARNING-- using %(dof)s degrees of freedom (%(dof)s dimensions)"%{"dof":self._degrees_of_freedom}
        slope /= (self._degrees_of_freedom*2.0)

        unit = "(%s)^2 / (%s)"%(self._distance_unit, self._time_unit)
        return(slope, unit)

    def find_diffusion_constant_in_units(self, begin, end, tounit="cm^2/s"):
        diffusion = self.find_diffusion_constant(begin, end)
        
        import os
        units_output = os.popen("units '%s' '%s' "%(str(diffusion[0])+" "+diffusion[1] , tounit)
                                , "r").read()
        #print units_output, "\n"
        new_diffcons = float(units_output.split("\n")[0].split()[1])
        return (new_diffcons, tounit)
