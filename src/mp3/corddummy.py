import logging ; thelog = logging.getLogger('mp3')
thelog.debug('Loading minimage.py')
import numarray
import mp3.cord
"""This sample of how to modify coordinates.

.
"""


class CordDummy(mp3.cord.Cord):
    """Minimal demonstratration of a cord-wrapper

    You can modify this to suit your needs. Note that error checking
    is on the sparse side.
    """
    def init(self):
        """Put any needed initilization here.

        Put your custom initialitilization stuff here. Do not remove
        this method, some other modules expect it to be here. Also,
        it must not break anything if it is called more than once.
        """
        framen = -1

    def setcord(self, cordobj):
        """Sets what we want to get our data from.

        Use this method to set the underlying coordinate-like object is
        (or whatever you need)
        """
        self.cord = cordobj
        self.cord.init()  # The cordobj can either be initilized or not, it should not break anything

    def nextframe(self):
        """
        """
        #          disp[0] -= (floorf((disp[0]/cstate->boxsize[0])+.5)/*-.5*/ ) * cstate->boxsize[0];
        frame = self.cord.nextframe()
        # do stuff to the frame
        self.frame = frame
        self.framen += 1
        return frame

    def zero_frame(self):
        """Returns to frame zero.
        """
        self.framen = -1
        self.cord.zero_frame()
        
    
    def read_n_frames(self, number):
        """Reads `number' more frames, only printing out the last one.

        """
        self.cord.read_n_frames(number-1)
        self.framen += number
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
        
