import logging ; thelog= logging.getLogger('mp3')
thelog.debug('Loading cords.py')
import mp3.cord

class CordMerge(mp3.cord.Cord):
    """A agregreator of coordinate sets.

    Internals: Call setcords() with a list of other cord objects.
    Then, reading from this new coordinate object will transparently
    concatenate them all together.
    """
    def __init__(self, cords=None):
        """
        """
        thelog.debug('--in cords.py, __init__()')
        # Editorial note: The new policy says that we always do all
        # initialization in the __init__ method, not by calling other
        # methods (unless there is a GREAT need).  So "cords" here
        # should never be none, really.
        if cords is not None:
            self._cords = cords
            self._init()
        else: 
            self._cords = [ ]
        
    def setcords(self, cords):
        """A wrapper to set the coordinates.

        Pass this method a list of other cord objects that you want this
        object to transparetly concatenate. This automatically calls the
        init() method.
        """
        thelog.debug('--in cords.py, setcords()')
        self._cords = cords
        self._init()

    def _init(self):
        """Sets up data needed for agregration.

        Calls init() on each of the underlying coordinate objects, then
        sets up needed variables like nframes and natoms.
        """
        thelog.debug('--in cords.py, init()')
        self._ncords = len(self._cords)
        #for cord in self.cords:
        #    cord.init()

        # We need some way to ensure that all of the internal dcds are
        # consistent.  Also, we need to store information about how many
        # frames each dcd has.

        # set up the number of atoms
        self._natoms = self._cords[0].natoms()
        for i in range(0,self._ncords):
            if self._cords[i].natoms() != self._natoms:
                thelog.error('natoms mismatch at dcd %d ', i)
        
        #set up the number of frames
        self._nframes = 0
        self._framelist = [ ]  # framelist[i] = cords[i].nframes
        self._cumulnframes = [ ] # index i has number of frames for all dcds up to cord i
        for i in range(0,self._ncords):
            self._framelist.append( self._cords[i].nframes())
            self._nframes += self._framelist[i]
            self._cumulnframes.append(self._nframes)

        # The following things are going to be set, without any real error checking (yet)
        # a normal dcd has these attributes:
        # ['framen', 'dcdtype', 'title', 'dcdfreq', 'tstep_size', 'frame', 'header_natoms', 'ntsteps', 'block_b', 'header_nframes', 'charm_v', 'initted', 'block_a', 'firsttstep', 'fo']
        

        #ntsteps should be dynamically generated...
            
        self._framen = -1
        self._initted = True
        
        


    def nextframe(self):
        thelog.debug('--in cords.py, nextframe()')
        self._framen += 1   #now framen records the frame that we need to spit back out
        thelog.debug('nextframe getting: %d', self._framen)
        for i in range(0, self._ncords):
            if self._framen < self._cumulnframes[i]:
                dcdnum = i
#                print self.cumulnframes[i]
#                print self.cords
                break
        self._frame = self._cords[dcdnum].nextframe()
        return self._frame

        
    def read_n_frames(self, number):
        """Reads this many frames.
        """
        thelog.debug("--cords.py, read_n_frames")
        if (number + self._framen) not in range(0, self._nframes):
            thelog.error("Supply a valid number for number of frames to read")
            thelog.error("  asked for %d more frames, but only %d total", number, self._nframes)
            thelog.error("")
            return None
        for i in range(0,number):
            self.nextframe()
        return self._frame

    def zero_frame(self):
        """Zero the frame location.

        Makes it where the first frame is the next one to be returned.
        (So you have to call nextframe() to get it)  This is framen = 0.
        """
        thelog.debug("--in cords.py, zero_frame()")
        thelog.debug("cords.py, zero_frame(), previous framen:%s"%self._framen)
        self._framen = -1
        for cord in self._cords:
            cord.zero_frame()
        
    
    def __getattr__(self, attrname):
        return self._cords[0].__dict__[attrname]
    
