#  This module provides an interface to aligning a coordinate
#  object as it is read
#
#
#

import logging ; thelog = logging.getLogger('mp3')
thelog.debug('Loading minimage.py')

import numarray
import Simplex
import mp3.functions
import mp3.cord


class CordAlign(mp3.cord.Cord):
    """Align coordinate objects as they are read.

    This class will, as it is reading coordinates, will translate/rotate
    the frame so that the msd of some atom subset is minimized with
    respect to the initial frame.
    
    Class contains:
        Methods:
            init() -- initlization, must be called manually after proper
                      setup done.
            setcord() -- sets where to get the coordinates to align
            setatoms() -- sets a possible atomlist
            nextframe_end_hook() -- hook right before nextframe() return. 
            

        Attributes:
            atomlist -- if != None, these atoms are used to do the alignment
            guess -- last guess of the minimum [x,y,z,  xangle,yangle,zangle]
            error -- error from the minimisation
            iterations -- number of iterations for last minimization
            firstframe -- (already sliced) frame we align to (exists only
                init() has been called)
            scale -- vector like "guess", but contains the initial scale
                of the simplex minimizer.
            
            _rmsd()  -- rmsd-wrapper that takes a guess vector and returns
                rmsd.  Handles slicing, too.
    """

    # properties of the transformation.  We save them for efficiency.
    guess = [0,0,0, 0,0,0]
    scale = (1,1,1, .1,.1,.1)

    # Atomlist tells which atoms to wrap by.  If it is None, then
    #  use all atoms.
    # nextframe_end_hook is called after each frame, for the user to
    #  replace with some other function to perform some arbitrary action.
    atomlist = None
    nextframe_end_hook = lambda x,y: True

    def __init__(self, cord, atoms=None):
        """Create the alignment object.

        The creation method can accept the following keyword arguments:
        cord= -- set where to get our coordinates
        atoms= -- If a list of atoms, use these atoms to align by.
                  If None, align by all atoms.

        If both of these keywords are given, the .init() method will be
        called for you.
        """
        # Set cord with the normal method if it was passed.  Set atoms if
        # atomlist was passed.  init() if both were passed.
        self.setcord(cord)
        self.setatoms(atoms)
        self.init()
    
    def init(self):
        """Initializes the aligner.

        This must be called manually after .setcord() and .setatoms() are
        used.
        """
        if self.atomlist is not None:
            self.firstframe = self.cord.nextframe()[self.atomlist].copy()
        else:
            self.firstframe = self.cord.nextframe().copy()
        self.cord.zero_frame()


        ### set up our state that depends on if we have an atomlist or not

        ## we either need to take a subset for the atoms, or we don't, for our rmsd function.
        ##  I'm going to implement it as an if right here.  I hypothesize that this will take less time
        ##  than having atomlist be all of it...

        ##  _rmsd a a function which wraps the RMSDing of the frame.
        #    self._rmsd = lambda vect, self: functions.rmsd(self.firstframe, \
        #             functions.cordtransform(self.curframe, move=vect[0:3], rotate=vect[3:6] ))
        #else:  # we don't have an atomlist
        #    self._rmsd = lambda vect, self: functions.rmsd(self.firstframe, \
        #             functions.cordtransform(self.frame, move=vect[0:3], rotate=vect[3:6] ))



    def setcord(self, cord):
        """Sets where to read our coordinates.

        You must call self.init() manually.
        """
        self.cord = cord

    def setatoms(self, atomlist=None):
        """Specifies which atoms we want to use to align our system

        Sets which atoms to align by.  If no argument is given or the argument
        is None, align by all atoms.

        You must call self.init() manually
        """
        self.atomlist = atomlist

        

    def nextframe(self):
        """Return the next frame after alignment
        """
        # extract the next frame, slice it by atoms if necessary.
        frame = self.cord.nextframe()
        if self.atomlist is not None:
            self.curframe = frame[self.atomlist].copy()
        else:
            self.curframe = frame.copy()
            
        # Create a wrapper function to take the msd between the two frames.
        # This is what is passed to the simplex optimizer.
        rmsd = lambda vect: mp3.functions.rmsd(self.firstframe, \
                    mp3.functions.cordtransform(self.curframe, move=vect[0:3], rotate=vect[3:6] ))
        self.minimizer = Simplex.Simplex(rmsd, self.guess, self.scale)
        self.guess, self.error, self.iterations = self.minimizer.minimize(monitor=0)
        self._frame = mp3.functions.cordtransform(frame, move=self.guess[0:3], rotate=self.guess[3:6])
        self.nextframe_end_hook(self)
        return self._frame


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

        # I will assume that any attributes not in self will be found in 
        # self.cord, and if not then I'll raise an exception.
        
        return self.cord.__dict__[attrname]
        
