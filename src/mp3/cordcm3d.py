import numarray
import mp3.cord

class CordCM3D(mp3.cord.Cord):
    """Reads CM3D coordinate trajectories.

    Creation:  CordCM3D(filename)

    It take a nontrivial amount of time to calculate the number of
    frames (we have to count each line, and there are a lot.)  There
    are methods to deal with this, documented below.  Other modules in
    mp3 won't handle this gracefully, so you want to figuer out how
    many frames there are total first, using either calcnframes() on
    setnframes(N).

    Methods:
        all inherited from Cord object.
        nextframe() -- like normal, raises ValueError if we run out of
            frames.
        nframes() -- returns None if we haven't calculated the number of
            frames yet.  Otherwise, behave as normal.
        calcnframes() -- Figure out the number of frames by counting
            up the number of lines in the file.  
        setnframes(N) -- If you want to tell the CordCM3D how many
            frames there are.

    The trajectory header looks like '# 15583 5 0.001'.  The first
    number is the number of atoms (natoms), the second one is stored
    in dcdfreq (how often a frame is saved), the last one is
    tstep_size (how long between timesteps).  In this example, a frame
    is stored every .005 picoseconds.
    """
    def __init__(self, confp):
        """Create the CordCM3D object.

        The keyword argument confp is the trajectory file filename.
        """
        self._filename = confp
        self._framen = -1
        fo = file(confp, "r")
        self._fo = fo
        line = fo.readline()
        #  '# 15583 5 0.001'
        line = line.split("#")[1]
        line = line.split()
        self._natoms = int(line[0])
        self._tstep_size = float(line[2])
        self._dcdfreq = int(line[1])
        self._frame = numarray.zeros(shape=(self._natoms,3), type=numarray.Float32)

    def nframes(self):
        """Return number of frames, or None if

        For CM3D, it is nontrivial to find the number of frames, so
        this may return None if the number of frames has not been
        calculated yet.
        """
        return getattr(self, "_nframes", None)

    def calcnframes(self, quick=True):
        """Count up the number of lines to find nframes
        """
        nlines = sum( [ 1 for line in file(self._filename) ] )
        # Quick method of calculating: calculate the number of lines.  
        nframes = (nlines-1) / (self._natoms+1)
        self._nframes = nframes
        #from rkddp.interact import interact; interact()
        return self._nframes
    def setnframes(self, nframes):
        """Set the total number of frames.

        This will simply save the number you give as the number of
        frames for this file.  Makes it easy if you already know it.
        """
        self._nframes = nframes
        

    def nextframe(self):
        """Loads the next frame into self.
        loads the next frame in self.frame . Hides details about wheter
        or not the dcd has been loaded into memory
        """
        nframes = self.nframes()
        # We must check if nframes is known yet.  If it isn't, then we
        # can't do this check.
        if nframes and self._framen >= nframes-1:
            log.critical("tried to read more frames than there were. BAILING OUT!!")
            return None
        frame = numarray.zeros(shape=(self._natoms,3), type=numarray.Float32)
        nextline = self._fo.readline  # function to directly call
        try:
            for i in range(self._natoms):
                # chokepoint
                frame[i,0], frame[i,1], frame[i,2] = [ float(X) for X in nextline().split() ]
            # read the unit cell parameters.
            # this may need to be refined some, if we want all parameters, not just orthogonal ones.
            cell = nextline()
            cell = cell.split()
            self._boxsize = (float(cell[0]), float(cell[4]), float(cell[8]))
        except ValueError:
            msg = "\n".join(["ValueError in reading the next frame.  This probably "
                             "means that we have ran out of frames in our file.  "+
                             "Last successful frame was %s, counting from zero.  "%self._framen+
                             "Current file is %s.  "%self._filename])
            mp3.log.error(msg)
            raise ValueError, msg
                           
        self._frame = frame
        self._framen += 1
        return self._frame


    def read_n_frames(self, ntoread):
        """Advance the current frame by the given number.

        Reads n frames from the CM3D trajectory file, only returning
        (and saving in self.frame) the last one.  If called with an
        argument of one, it behaves exactly like nextframe().  This can't go backwards.
        """
        if ntoread < 1:
            raise Exception, "We must read at least one frame."
        nextline = self._fo.readline
        try:
            for i in range(ntoread-1):
                for j in range(self._natoms):
                    nextline()
                nextline()
                self._framen += 1
        except ValueError:
            msg = "\n".join(["ValueError in reading the next frame.  This probably "
                             "means that we have ran out of frames in our file.  "+
                             "Last successful frame was %s, counting from zero.  "%self._framen+
                             "Current file is %s.  "%self._filename])
            mp3.log.error(msg)
            raise ValueError, msg
        return self.nextframe()

    def zero_frame(self):
        """makes it where the next frame returned is 0. This function itself
        does not return anything.
        """
        del self._fo
        # we'll just re-initilize it all.  That can't hurt.  And it
        # doesn't slow it down too much.
        self.__init__(self._filename)
