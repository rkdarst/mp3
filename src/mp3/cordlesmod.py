import logging ; thelog = logging.getLogger('mp3')
thelog.debug('Loading minimage.py')
import numarray
import mp3.cord
import mp3.functions
"""This does un-wrapping of LES'ed cord objects.

LES = locally enhanced sampling


"""


class CordLESMod(mp3.cord.Cord):
    """Minimal demonstratration of a cord-wrapper

    You can modify this to suit your needs. Note that error checking
    is on the sparse side.
    """
    #def init(self):
    #    """Put any needed initilization here.

    #    """
    #    pass

    #
    #  "Initilization"
    #

    def __init__(self, **keywords):
        """
        """
        self._which_frame = -1      # by default, return the average frame 
        self._natoms = -1           # temporary identifier saying that natoms
                                    #  hasn't been calculated yet.

        if keywords.has_key("cord"):
            self.setcord(keywords["cord"])
        if keywords.has_key("numrep"):
            self.setnumrep(keywords["numrep"])
        if keywords.has_key("whichframe"):
            self.setwhichframe(keywords["whichframe"])
        if keywords.has_key("leslist"):
            self.setleslist(keywords["leslist"])
        if keywords.has_key("atomlist"):
            self.setatomlist(keywords["atomlist"])


    def setcord(self, cord):
        """Sets what we want to get our data from.
        Use this method to set the underlying coordinate-like object is
        (or whatever you need)
        """
        self.cord = cord


    def setnumrep(self, numrep):
        """Sets the number of replicates for each replicated atom.
        """
        self._numreplicates = numrep


    def setatomlist(self, atomlist):
        # This must be called last
        self._atomlist = atomlist
        print "calculating leslist"
        self._leslist = mp3.functions._atomlist_to_leslist(atomlist, self.natoms())
        print "done calc leslist"


    def setleslist(self, leslist):
        """Set the list describing the locally enhanced sampling in use.
        """
        self._leslist = leslist


    def setwhichframe(self, whichframe):
        """Specify what self.frame() should be.

        This module finds both the average frame and the frame of each
        replicate for the LES.

        -1   ==> the average
        >= 0 ==> which replicate's frame to return.
        """
        self._which_frame = whichframe

    #
    #   "attributes"
    #


    def frame(self):
        return self._frame

    def frames(self):
        """Matrix containing all frames:

        self.frames()[which_replicate, atom_number, coordinate]
           (of course, numbering starts from zero)
        """
        return self._frames

    def mean(self):
        if not hasattr(self, "_mean"):
            self._calculate_mean()
        return self._mean

    def variance(self):
        if not hasattr(self, "_variance"):
            self._calculate_variance()
        return self._variance

    def rmsd(self):
        if not hasattr(self, "_rmsd"):
            self._calculate_rmsd()
        return self._rmsd

    def natoms(self):
        if self._natoms == -1:
            self._calculate_natoms()
        return self._natoms

    def _calculate_natoms(self):
        """Figure out how many atoms our final system has.

        This is needed when you specify an atomlist (list each
        replicate) and a cord object (has the number of atoms,
        including replicates).  We need to figure out how many
        atoms there are in the system from this.

        1) If self._leslist is specified, then natoms = sum(leslist)
          - must have self._leslist
        If we don't have leslist, then we need another plan.
        In order to make leslist, we must first have natoms.  So we
        _must_ have an atomlist in that case..

        2) Otherwise, use
        natoms = cord.natoms() - len(atomlist)*(numrep-1)
          - must have self.cord
          - must have self._atomlist
          - must have self._numreplicates
        """
        if hasattr(self, "_leslist"):
            self._natoms = sum(self._leslist)
            return self._natoms
        elif hasattr(self, "cord") \
                and hasattr(self, "_atomlist") \
                and hasattr(self, "_numreplicates"):
            self._natoms = self.cord.natoms() \
                           - (len(self._atomlist) * (self._numreplicates-1))
            return self._natoms
        else:
            print "Error !"
            mp3.log.error("Can not calculate natoms with given information!")


    def nextframe(self):
        """Return the next frame.
        """

        self.cord.nextframe()
        self._clear_cache()
        
        # Some convenience variables
        leslist = self._leslist
        cord = self.cord
        numreplicates = self._numreplicates
        
        frame = numarray.zeros((self.natoms(), 3), type=numarray.Float32)
        oldframe = cord.frame()
        position = 0  # spot in leslist
        newpos = 0    
        oldpos = 0

        # Iterate through leslist, using "break" to escape.
        while True: 
            if len(leslist) > position:
                length = leslist[position]      # lenght is how long this block is.
                #print frame[newpos:newpos+length].shape
                #print oldframe[oldpos:oldpos+length].shape
                frame[newpos:newpos+length] = oldframe[oldpos:oldpos+length]
                newpos += length
                oldpos += length
                position += 1
            else:
                break
            if len(leslist) > position:
                length = leslist[position]
                newpos += length 
                oldpos += length * numreplicates
                position += 1
            else:
                break

        # So, at this point we have filled "frame" with just the non-
        # replicated properties.  Now, we have to go back through and
        # fill it with the atoms that have been replicated.

        frames = numarray.zeros((numreplicates, self.natoms(), 3),
                                type=numarray.Float32)

        for i in range(numreplicates):
            frames[i] = frame.copy()
        position = 0
        newpos = 0
        oldpos = 0
        while True:
            if len(leslist) > position:
                length = leslist[position]
                newpos += length
                oldpos += length
                position += 1
            else:
                break
            if len(leslist) > position:
                length = leslist[position]
                for k in range(length):
                    for l in range(numreplicates):
                        frames[l,newpos] = cord.frame()[oldpos].copy()
                        oldpos += 1
                    newpos +=1
                position +=1
            else:
                break
        #print "we got past making the arrays"
        self._frames = frames

        # Now, pick which one to return:
        if self._which_frame == -1:
            self._frame = self.mean()
        else:
            self._frame = self._frames[self._which_frame]
        return self._frame

    def _clear_cache(self):
        if hasattr(self, "_mean"):
            del self._mean
        if hasattr(self, "_variance"):
            del self._variance
        if hasattr(self, "_rmsd"):
            del self._rmsd


    def _calculate_mean(self):
        # Now that we have filled frames with each LES thing, we can
        # proceed to calculate the averages and all.
        mean = numarray.sum(self._frames, axis=0) / self._numreplicates
        self._mean = mean
    def _calculate_variance(self):
        variance = numarray.subtract(self._frames, self.mean())
        variance *= variance
        variance = numarray.sum(variance, axis=0)
        variance = numarray.sqrt(variance)
        variance /= self._numreplicates
        self._variance = variance
    def _calculate_rmsd(self):
        rmsd = numarray.subtract(self._frames, self.mean())
        rmsd *= rmsd
        rmsd = numarray.sum(rmsd, axis=2)
        rmsd = numarray.sum(rmsd, axis=0)
        rmsd /= self._numreplicates
        numarray.sqrt(rmsd, rmsd)
        self._rmsd = rmsd
    


    def zero_frame(self):
        """Returns to frame zero.
        """
        self.cord.zero_frame()
        
    
    def read_n_frames(self, number):
        """Reads `number' more frames, only printing out the last one.

        """
        self.cord.read_n_frames(number-1)
        return self.nextframe()


    def __getattr__(self, attrname):
        """Wrapper for getting cord attributes.

        This is a wrapper for getting things like self.nframes when you
        don't have to worry about setting them yourself.  Going on the
        hypothesis that most of the time these aren't changed, for
        anything that you haven't defined yourself, it will pass it
        through to self.cord.

        Note that this could be bad in some cases!  But I'll take care
        of them when I find them.
        """
        return self.cord.__dict__[attrname]
        
