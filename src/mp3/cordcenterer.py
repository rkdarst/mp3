import logging ; thelog = logging.getLogger('mp3')
thelog.debug('Loading minimage.py')
import numarray
import mp3.cord
"""This module does basic min

.
"""


class CordCenterer(mp3.cord.Cord):
    """
    """

    def __init__(self, **keywords):
        """Create the centerer object.

        The creation method can accept the following keyword arguments:
        cord= -- set where to get our coordinates.
        weights= -- Should be a list of all atoms' weights to align by.

        If both of these keywords are given, the .init() method will be
        called for you.  Both of these keywords are optional, if they
        are not given, you must use the .setcord() and .setweights()
        methods to set this information.
        """
        # Set cord with the normal method if it was passed.  Set atoms if
        # atomlist was passed.  init() if both were passed.
        if keywords.has_key("cord"):
            self.setcord(keywords["cord"])
        if keywords.has_key("weights"):
            self.setweights(keywords["weights"])
        if keywords.has_key("cord") and keywords.has_key("weights"):
            self.init()

    def init(self):
        """Initilize the centerer."""
        self.cord.init()
        
    def setcord(self, cord):
        """Sets what we want to min image wrap.

        Use this method to set the underlying coordinate-like object is
        (or whatever you need)
        """
        self.cord = cord
        #self.cord.init()  # The cordobj can either be initilized or not, it should not break anything

    def setweights(self, weights):
        """Tells it which atoms you want center.

        Pass this method a (n,) numarray listing the relative weights of the atoms,
        and this will translate the thing such that the center of mass of these
        atoms is at the origin.
        """
        # I call it weights because it is the relative weights of the atoms, not
        #  because I'm confusing weight and mass
        self.weights = weights
        self.totalweight = numarray.sum(weights)
        

    def nextframe(self):
        """
        """
        #          disp[0] -= (floorf((disp[0]/cstate->boxsize[0])+.5)/*-.5*/ ) * cstate->boxsize[0];
        frame = self.cord.nextframe().copy() #don't want to mess up the trasposing
        frame.transpose()
        newframe = numarray.multiply(frame, self.weights)
        sums = numarray.sum(newframe, axis=1)
        numarray.divide(sums, self.totalweight, sums)
        frame = numarray.subtract(self.cord.frame, sums)
        self._frame = frame
        return self._frame

    def zero_frame(self):
        self.cord.zero_frame()
    
    def read_n_frames(self, number):
        self.cord.read_n_frames(number-1)
        return self.nextframe()

    def __getattr__(self, attrname):
        return self.cord.__dict__[attrname]

        
