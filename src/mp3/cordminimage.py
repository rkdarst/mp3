import logging ; thelog = logging.getLogger('mp3')
thelog.debug('Loading minimage.py')
import numarray
import mp3.cord
#import mp3.xst
"""This module does basic minimum image wrapping."""


class CordMinImage(mp3.cord.Cord):
    """Do minimum image wrapping.

    This object will translate coordinates by wrapping all values to
    fit inside of some given box.  
    """
    def init(self):
        """Not needed.

        Exists so that calling it won't fail.
        """
        pass

    def __init__(self, **keywords):
        """Initilize the coordinate object.

        You can use the keyword arguments "cord" and "boxsize" to set
        the coordinates and boxsize.  For the format of the boxsize
        argument, see the setboxsize() method.

        If both "cord" and "boxsize" are given here, the method init()
        will be called.
        """
        if keywords.has_key("cord"):
            self.setcord(keywords["cord"])
        if keywords.has_key("boxsize"):
            self.setboxsize(keywords["boxsize"])
        if keywords.has_key("cord") and keywords.has_key("boxsize"):
            self.init()
        
    def setboxsize(self, boxsize):
        """Sets what boxsize you want to do minimum image wraping for.

        Argument should be a tuple, list, or numarray of the form [x,y,z].
        """
        # Here is how we will do this.  self._boxsize_call will be a 
        # callable which will return the new boxsize.  Thus, this could 
        # be an xst object which when called will return the next boxsize.
        # If we have a constant box size, we will use a lambda function
        # which always returns the same boxsize.

        # This method will accept either a tuple or callable (such as an
        # xst object).  If the argument is callable, we just use it as
        # the function to return the next boxsize.  Otherwise, we 
        # construct a labmda function to always return the argument.

        if callable(boxsize) is True:
            self._boxsize_call = boxsize
        else:
            if len(boxsize) != 3:
                thelog.error("Given boxsize does not have three dimensions!")
                raise   # FIXME
            boxsize_call = lambda : boxsize    # this is broken, we need a lambda function which
                                                 # takes no arguments
            self._boxsize_call = boxsize_call

    def setcord(self, cordobj):
        """Sets what we want to min image wrap.

        Use this method to set the underlying coordinate object.
        """
        self.cord = cordobj

    def nextframe(self):
        """Returns the next frame.

        This will return the next frame, wrapping it.
        """
        self._frame = self.cord.nextframe()
        self._wrap_frame

    def _wrap_frame(self):
        self.boxsize = numarray.asarray(self._boxsize_call())
        frame = self._frame
        wrapping = frame.copy()
        numarray.divide(wrapping, self.boxsize, wrapping)
        numarray.add(wrapping, .5, wrapping)
        numarray.floor(wrapping, wrapping)
        numarray.multiply(wrapping, self.boxsize, wrapping)
        numarray.subtract(frame, wrapping, frame)
        self._frame = frame
        return _frame
            
    def nextframe_nowrap(self):
        """Advance the frame without wrapping."""
        self._frame = self.cord.nextframe()

    def read_n_frames(self, number):
        self.cord.read_n_frames(number-1)
        return self.nextframe()

    def __getattr__(self, attrname):
        if attrname not in self.__dict__.keys():
            return self.cord.__dict__[attrname]
        else:
            return self.__dict__[attrname]
        
