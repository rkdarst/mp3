#Use only a subset of atoms in a DCD object.
#
#This module defines a class that will extract a (constant) subset of
#atoms from a DCD object, and only pass those atoms on.
#

import logging ; thelog = logging.getLogger('mp3')
thelog.debug('Loading cordatomslice.py')
import numarray
import mp3.cord


class CordAtomSlice(mp3.cord.Cord):
    """Use only a subset of atoms in a cord object.
    
    This module defines a class that will extract a (constant) subset of
    atoms from a cord object, and only pass those atoms on.  This might
    be useful when, for instance, you want to process a DCD of just a
    protein instead of the protein + solvent.

    You can initilize the instance with the optional keyword arguments
    'cord' and 'atoms'.  These will set the cord object to get
    coordinates from, and which atoms to pass through.  If both are
    given, then .init() method will be called for you.
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
        if keywords.has_key("cord"):
            self.setcord(keywords["cord"])
        if keywords.has_key("atoms"):
            self.setatoms(keywords["atoms"])
        if keywords.has_key("cord") and keywords.has_key("atoms"):
            self.init()

    
    def init(self):
        """Initilize the atomslicer."""
        self.cord.init()

    def setcord(self, cordobj):
        """Set the coordinate object to get our data from.

        Use this method to set the underlying coordinate-like object is.
        The slice is taken from this coordinate object.
        """
        self.cord = cordobj
        self.cord.init()  # The cordobj can either be initilized or not, it should not break anything

    def setatoms(self, atomlist):
        """Set which atoms to pass through.

        The argument to this function is a list of atoms which will be
        passed through the slicer.  There is no default to use all
        atoms, since if that was likely, then this wrapper wouldn't be
        used, would it ?
        """
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
        
