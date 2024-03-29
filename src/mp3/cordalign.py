#  This module provides an interface to aligning a coordinate
#  object as it is read
#
#
#
import sys
import logging ; thelog = logging.getLogger('mp3')
thelog.debug('Loading minimage.py')

import numpy
import mp3.log
try:
    import scipy.optimize
except:
    mp3.log.debug("We did not load scipy.optimize library, CordAlign will be broken.")
    pass
import mp3.functions
import mp3.cord


class CordAlign(mp3.cord.Cord):
    """Align coordinate objects as they are read.

    This class will, as it is reading coordinates, will translate/rotate
    the frame so that the msd of some atom subset is minimized with
    respect to the initial frame.
    
    Class contains:
        Methods:
            __init__() -- See docstring for initilization
                          information.
            nextframe() -- Do alignment and return the next aligned frame.
                          Aligned frame also stored in frame() method.
            nextframe_end_hook() -- hook right before nextframe() return.
            averageframe() -- return average of atom positions after
                          alignment (see docstring)

        Attributes:
            atomlist -- if != None, these atoms are used to do the alignment
            minimizer -- the minimization algorithm for the alignment
            guess -- last guess of the minimum [x,y,z,  xangle,yangle,zangle]
            error -- error from the minimisation
            iterations -- number of iterations for last minimization
            verbose -- increse verbosity if true
            aligntoframe -- Frame we align to.  Only contains the atoms we
                     are using for aligning.
            scale -- vector like "guess", but contains the initial scale
                of the simplex minimizer.
            _rmsd()  -- rmsd-wrapper that takes a guess vector and returns
                rmsd.  Handles slicing, too.

    Note: All atomlists should be sorted before using them here.  This
    is true in general with mp3.
    """

    # properties of the transformation.  We save them for efficiency.
    # We start off with these default values.
    guess = [0,0,0, 0,0,0]
    scale = (1,1,1, .1,.1,.1)

    # Atomlist tells which atoms to wrap by.  If it is None, then
    #  use all atoms.
    # nextframe_end_hook is called after each frame, for the user to
    #  replace with some other function to perform some arbitrary action.
    nextframe_end_hook = lambda x,y: True

    def __init__(self,
                 cord,
                 atomlist=None,
                 aligntoframe=None,
                 saveaverage=False,
                 minimizer="scipy:powell",
                 callback=None,
                 verbose=False):
        """Create the alignment object.

        The creation method can accept the following keyword arguments:
        cord=
            Our child Cord object, we get the coordinates to align from
            here.
        atomlist=
            If a list of atoms, use these atoms to align by.  If None,
            align by all atoms.
        aligntoframes=
            If this is None, then align future frames to the first frame
            (using the given atomlist).  Otherwise, this is a frame of
            _just_ the atoms to align by (ex: aligntoframe=
            whole_frame[atomlist] ).  You should specify an 'atomlist'
            consistent with this frame (same number of atoms, or None for
            the whole frame).  
        saveaverage=
            If we specify True here (default False), then we will save
            the data needed to find an average frame at some point.  See
            the docstring the averageframe() method.
        minimizer=
            Specifies the minimization routine to be used in aligning the
            structure.  The default is the direction set (powell's) method.
            Another option is the downhill simplex (Nelder-Mead) method.
        callback=
            A function which will be called on each simplex minimization
            step.  It is called as `callback(xk, self)` with `xk` the
            current guess and `self` being the CordAlign object.
        verbose=
            The default state is quiet.  With verbose specified as True, information
            from minimization routine will be output.  This includes the function value,
            the number of iterations, and an error estimate.
        """
        # Set cord with the normal method if it was passed.  Set atoms
        # if atomlist was passed.  (These came from the old
        # .setatoms() and setcord() methods)
        self.cord = cord
        self.atomlist = atomlist
        self._align_to_next_frame = False
        self.minimizer = minimizer
        self.callback = callback
        self.verbose = verbose

        self._saveaverage = saveaverage
        if saveaverage:
            self._sum_of_frames = numpy.zeros(dtype=numpy.float32,
                                                 shape=(self.cord.natoms(), 3))
        self._sum_of_frames_count = 0
        
        # Some of this stuff came from the old .init() method.
        if aligntoframe is not None:
            # if we are given aligntoframe:
            # test to be sure that we have the right number of atoms.
            if self.atomlist is not None and aligntoframe.shape[0] != len(self.atomlist):
                mp3.log.error("You specified a aligntoframe to and an atomlist.")
                mp3.log.error(" But the frame didn't have the same number of atoms")
                mp3.log.error(" as the atomlist, thus I don't know what to do here.")
                mp3.log.error("The aligntoframe should already be reduced to only")
                mp3.log.error(" include the atoms that are already in the atomlist.")
                mp3.log.error("The error will appear later in the code.")
            if self.atomlist is None and aligntoframe.shape[0] != self.cord.natoms():
                mp3.log.error("You specified a aligntoframe without an atomlist,")
                mp3.log.error(" so we default to all atoms.  But the number of atoms")
                mp3.log.error(" in aligntoframe isn't the same as the number of atoms")
                mp3.log.error(" in the cord.  This is wrong.")
                mp3.log.error("The error will appear later in the code.")
            self.aligntoframe = aligntoframe
            # self.atomlist should already be set correctly, as per above.
        else:
            # if we are not given aligntoframe, so we want to align by the first frame:
            self._align_to_next_frame = True
            # We store it here, and then act on it on the first read.
            # We don't want to set self.aligntoframe now, because that
            # would require us to read forward, then back up again,
            # which is isn't perfectly supported right now, and
            # probably isn't the best thing to be doing anyway.

        ### set up our state that depends on if we have an atomlist or not

        ## we either need to take a subset for the atoms, or we don't,
        ##  for our rmsd function.  I'm going to implement it as an if
        ##  right here.  I hypothesize that this will take less time
        ##  than having atomlist be all of it...

        ##  _rmsd a a function which wraps the RMSDing of the frame.
        #    self._rmsd = lambda vect, self: functions.rmsd(self.firstframe, \
        #             functions.cordtransform(self.curframe, move=vect[0:3], rotate=vect[3:6] ))
        #else:  # we don't have an atomlist
        #    self._rmsd = lambda vect, self: functions.rmsd(self.firstframe, \
        #             functions.cordtransform(self.frame, move=vect[0:3], rotate=vect[3:6] ))

    def nextframe(self):
        """Return the next frame after alignment

        Uses algorithm and notation from Kabsch, Wolfgang, (1976) "A solution of the best 
        rotation to relate two sets of vectors", Acta Crystallographica 32:922
        
        """
        # extract the next frame, slice it by atoms if necessary.
        frame = self.cord.nextframe()
        self.ca_origframe = frame  # not .copy()'ing it, when I wrote it we didn't need to.
        if self.atomlist is None:
            self.curframe = frame.copy()
        else:
            self.curframe = frame[self.atomlist].copy()

        if self._align_to_next_frame:
            # We wanted to align to the first frame of the DCD.  See
            # above for an explanation.
            if self.atomlist is None:
                self.aligntoframe = frame.copy()
            else:
                self.aligntoframe = frame[self.atomlist].copy()
            self._align_to_next_frame = False

        X = self.aligntoframe.copy()
        Y = self.curframe.copy()

        natoms,ndimensions = numpy.shape(X)
        
        center1 = sum(X,0) / float(natoms)
        center2 = sum(Y,0) / float(natoms)
        X -= center1
        Y -= center2

        E0 = sum(sum(X * X)) + sum(sum(Y * Y))

        correlation_matrix = numpy.dot(numpy.transpose(Y), X)

        V, S, W_trans = numpy.linalg.svd(correlation_matrix)

        is_reflection = (numpy.linalg.det(V) * numpy.linalg.det(W_trans)) < 0.0
        if is_reflection:
            # reflect along smallest principal axis
            S[-1] = -S[-1]
            V[-1,:] = V[-1,:] * (-1.0)

        optimal_rotation = numpy.dot(V, W_trans)
        self._frame = numpy.dot(frame, optimal_rotation) - center2 + center1
        
        self.nextframe_end_hook(self)
        return self._frame
              
        
        # UPDATE (JRD, Sept. 2007):  This section of code contained between the #**...*# markers
        # is not the correct way to align frames for a single protein MD trajectory.  The proper 
        # method for aligning frames is due to Kabsch:  Kabsch, Wolfgang, (1976) "A solution of 
        # the best rotation to relate two sets of vectors", Acta Crystallographica 32:922
        #*******************************************************************************************#
        # Create a wrapper function to take the msd between the two frames.
        # This is what is passed to the simplex optimizer.
        #rmsd = lambda vect: mp3.functions.rmsd(self.aligntoframe, \
        #            mp3.functions.cordtransform(self.curframe, move=vect[0:3], rotate=vect[3:6] ))
        #if self.verbose:
        #    dispval = 1
        #else:
        #    dispval = 0
        #if self.callback:
        #    callback = lambda xk:  self.callback(xk, self)
        #    # your callback can increment this to figure out what step it is on
        #    self.iterstep = 0 
        #else:
        #    callback = None
        #if self.minimizer == "scipy:powell":
        #    result = scipy.optimize.fmin_powell(rmsd,self.guess,disp=dispval,full_output=1,
        #                                        ftol=1e-6,callback=callback)
        #    self.iterations = result[3]
        #    self.funcalls = result[4]
        #elif self.minimizer == "scipy:simplex":
        #    result = scipy.optimize.fmin(rmsd,self.guess,disp=dispval,full_output=1,
        #                                 callback=callback)
        #    self.iterations = result[2]
        #    self.funcalls = result[3]
        #else:
        #    sys.stderr.write("ERROR: minimizer must be either scipy:powell or scipy:simplex")
        #    sys.exit()
        #self.guess = result[0]
        #self._frame = mp3.functions.cordtransform(frame, move=self.guess[0:3],
        #                                          rotate=self.guess[3:6])
        #if self._saveaverage:
        #    self._sum_of_frames += self._frame
        #    self._sum_of_frames_count += 1
        #
        #self.nextframe_end_hook(self)
        #return self._frame
        #*******************************************************************************************#        

    def zero_frame(self):
        """Returns to frame zero.
        """
        self.cord.zero_frame()
    
    def read_n_frames(self, number):
        """Reads `number' more frames, only printing out the last one.

        """
        self.cord.read_n_frames(number-1)
        return self.nextframe()

    def averageframe(self):
        """Returns the average of all so-far aligned frames.

        As we go through the trajectory, we save the data needed to
        compute the average of the whole frame.  This method returns
        the average of all atoms' locations, in a simple numpy array
        (looks just like a normal frame).  The average is taken for
        each atom, not just those being used to align.  Some averages
        may not be relavent, for example if you are aligining a
        protein in water by it's backbone (protein averages are
        meaningful, but not for water since they are randomly
        diffusing around).

        In order to use this, you must have initilized the CordAlign
        with saveaverage=True.
        """
        if self._saveaverage:
            averageframe = self._sum_of_frames / self._sum_of_frames_count
            return averageframe
        else:
            mp3.log.error("Cannot take average frame.  Object initilized")
            mp3.log.error(" with saveaverage=False, so we haven't been saving data.")

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
        
