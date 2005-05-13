#Use only a subset of atoms in a DCD object.
#
#This module defines a class that will extract a (constant) subset of
#atoms from a DCD object, and only pass those atoms on.
#

import numarray
import mp3.cord

from mp3.log import mp3log
mp3log.debug("loading cordatomslice.py")

class CordAtomSlice(mp3.cord.Cord):
    """Use only a subset of atoms in a cord object.

    Example usage:  (protein with solvent) -> (just the protein)

    >>> C = mp3.CordDCD("my.dcd")
    >>> C.nframes()
    1000
    >>> C.natoms()
    100000
    
    >>> C_sliced = mp3.CordAtomSlice(cord=C, atomlist=range(3000))
                                   # the first 3000 atoms are the protein
    >>> C_sliced.nframes()
    1000
    >>> C_sliced.natoms()
    3000

    
    This module defines a class that will extract a (constant) subset
    of atoms from a cord object, and only pass those atoms on for any
    further analysis or processing.  This might be useful when, for
    instance, you want to analyze the coordinates of just a protein
    instead of the protein + solvent.

    The class should be initilized as such:

    C = CordAtomSlice(cord=CORD, atomlist=ATOMLIST)

    CORD should be some other cord object.

    ATOMLIST is a specification of which atoms to include.
    Technically, the actual slicing is performed this way:

            frame = nextframe[ATOMLIST]

    So ATOMLIST is anything that can "slice" a numarray.  Most typically,
    this would be a simple atomlist:

    atomlist = range(10)     # the first 10 atoms

    atomlist = [9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
                             # the first 10 atoms IN REVERSE ORDER.
                             # This is a feature of numarray slicing.
                             # It will return the atoms in the order
                             # you specify.

    atomlist = [0, 1, 8, 9]
                             # Atoms 1,2, 9, 10 (note counting from zero)

    atomlist = list(S.findatoms( <stuff> ))
    atomlist.sort()
                             # If you use findatoms, remember to
                             # convert it to a list first.  Also, you
                             # probably want to sort it, or else the
                             # indexes may be scrambled and you can't
                             # use it

    """

    def __init__(self, **keywords):
        """Create the atomslice object.

        The creation method can accept the following keyword arguments:
        cord= -- set where to get our coordinates.
        atoms= -- Should be a list of atoms to align by.

        If both of these keywords are given, the .init() method will be
        called for you.  Both of these keywords are optional, if they
        are not given, you must use the .setcord() and .setatoms()
        methods to set this information.
        """
        # Set cord with the normal method if it was passed.  Set atoms if
        # atomlist was passed.  init() if both were passed.
        mp3log.debug("Initilizing CordAtomSlice instance.")

        keywords.has_key("atomlist") or mp3log.warn(
            "Not initializing CordAtomSlice with an atomlist.")
        keywords.has_key("cord") or mp3log.warn(
            "Not initializing CordAtomSlice with an cord.")
        
        if keywords.has_key("cord"):
            self.setcord(keywords["cord"])

        if keywords.has_key("atomlist"):
            self.setatomlist(keywords["atomlist"])

        if keywords.has_key("cord") and keywords.has_key("atoms"):
            #self.init()
            pass

    
    def init(self):
        """Initilize the atomslicer."""
        pass
        #self.cord.init()

    def setcord(self, cord):
        """Set the coordinate object to get our data from.

        Use this method to set the underlying coordinate-like object is.
        The slice is taken from this coordinate object.
        """
        self.cord = cord

        # We don't use the init() method anymore, the instances should
        # be initilized as soon as they are created.
        #self.cord.init()  
                         

    def setatomlist(self, atomlist):
        """Set which atoms to pass through.

        The argument to this function is a list of atoms which will be
        passed through the slicer.  There is no default to use all
        atoms, since if that was likely, then this wrapper wouldn't be
        used, would it ?
        """
        mp3log.debug("We have %s atoms passed through."%len(atomlist))
        self.atomlist = atomlist
        self._natoms = len(atomlist)

    def nextframe(self):
        self._frame = self.cord.nextframe()[self.atomlist].copy()
        return self._frame

    def zero_frame(self):
        self.cord.zero_frame()
        
    def read_n_frames(self, number):
        self.cord.read_n_frames(number-1)
        return self.nextframe()

    def __getattr__(self, attrname):
        return self.cord.__dict__[attrname]
        
