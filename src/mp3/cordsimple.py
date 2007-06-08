# Richard Darst, 2005
#
#

import mp3.cord


class CordSimple(mp3.cord.Cord):
    """Simple way to make Cord objects.

    Most Cord objects read their data from files.  This works well
    except for the fact that you can't modify the coordinates very
    well, or create your own coordinates.  CordSimple gets around
    these limitations with a very simple interface.

    All data is stored in a list (CordSimple.all_frames).  Each
    element of the list is a numpy array representing one frame.  Element
    zero is frame zero, etc.  You may manipulate this list hovever you
    wish.  Note that there is not much error checking, so it won't
    complain if your frames in the list don't have the same number of
    atoms.  Instead, it will break later and you will be left
    scratching your head.

    Methods such as CordSimple.nframes() or CordSimple.natoms() are
    found by examining this list, for example these are
    len(CordSimple.all_frames) for nframes() and
    len(CordSimple.all_frames[0]) for natoms().

    The framen() and nextframe() ideas are still valid.  To handle
    this, we maintain an internal counter telling which is our current
    frame.  Asking for a frame returns that element from the list.
    """

    def __init__(self):
        self._framen = -1 
        self.all_frames = [ ]
    # Our primary data structure here is self._all_frames.  It is as
    # list, each element of the list is a frame.  We get as much
    # information from this list and the frames in it as we possibly
    # can.  That means we generate it on the fly.
    #
    def nframes(self):
        """Number of frames.

        return len(self.all_frames)
        """
        return len(self.all_frames)
    def natoms(self):
        """Number of atoms.

        return len(self.all_frames[0])
        """
        return len(self.all_frames[0])
    def frame(self):
        """Return the current frame.
        
        return self.all_frames[self.framen()]
        """
        return self.all_frames[self.framen()]
    def nextframe(self):
        """Return next frame.
        
        self._framen += 1
        return self.frame()
        """
        self._framen += 1
        return self.frame()

    def read_n_frames(self, n):
        self._framen += 1
        return self.frame()
    def zero_frame(self):
        self._framen = -1
    
    def appendframe(self, frame):
        """Append the given frame onto our frame list.
        
        self.all_frames.append(frame)
        """
        self.all_frames.append(frame)
