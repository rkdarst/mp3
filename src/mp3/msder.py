

import numpy
import mp3.functions

class Msder:
    """Perform mean-square calculations on coordinate data.

    What is calculated is the distribution of
       mean( (atomA_frame1 - atomA_frame2)**2 )
    for all pairs of atomA's in 'atomlist', as a function of the
    distance between frame1 and frame2.

    This class is initilized with the keyword arguments 'cord', 'window',
    and optionally 'atomlist'
      cord == The coordinate object to get data from.  We begin gathering
              data from the next frame in the object.
      window == The maximum distance between frames to be sampled.  Units
              are frames.
      atomlist == list contaning which atoms to analyze.  If not given or
              None, all atoms will be used        


    Analysis is done 'by frame'.  When you first initilize the object,
    no analyzes are done.  Each time you call the do_n_msd(n) method,
    it will compare the first frame to each of the next 'window'
    frames, tabulating data.
    
    Note that we get an equal amount of data for all distances, as in
    we don't do more calculations for the short distances.  If this is
    you goal, use the 'finish_msd_window' method after you have read
    as many frames as you desire.
    
    The 'msd' method returns a numpy array with the msd information.
      msdlist[i] == separation of i frames apart.
    Note that this means that msdlist[0] == 0, since being zero frames
    apart mean taking the deviation of a frame with itself.
        
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

        You can use the method 'do_all_msd' to automatically do a MSD
        calculation for an entire cord object.  It assumes that you
        hand it a fresh cord object, with no frames read yet (it will
        do them 'all'.)

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
        """Return a numpy array containing the msd.
        """
        return numpy.divide(self._msd_sum, self._count_sum)

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
