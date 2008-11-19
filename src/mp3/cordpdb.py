# This module serves as an interface to PDBs.  It allows you to
# access the coordinates stored in PDBs.
#
#

import logging ; thelog = logging.getLogger('mp3')
thelog.debug('Loading pdbcord.py')
import glob
import numpy
import mp3.cord


class CordPDB(mp3.cord.Cord):
    """Manages coordinates from a PDB series
    """
    initted = False

    def __init__(self, pdblist=None):
        """
        """
        if pdblist != None:
            self.setpdblist(pdblist)

    def init(self):
        """Initilizes data about the PDBs.

        You don't need to use this, use self.setpdblist and this will
        be called automatically.
        """
        # check to see if we have already been initilized, if so, don't
        # do it again.
        if self._initted == True:
            return None
        else:
            self._initted = True

        # find number of frames and number of atoms
        self._nframes = len(self._pdblist)
        self._framen = -1
        # number of atoms
        # we're going to count up how many lines we have
        # that start with "ATOM".  Once we get to a line
        # that doesn't start with that, we stop.
        # Leading garbage is discarded.
        self._natoms = 0
        pdbfo = file(self._pdblist[0], "r")
        while True:
            line = pdbfo.readline()
            if line[0:4] == "ATOM":
                break
        self._natoms += 1
        while True:
            line = pdbfo.readline()
            if line[0:4] == "ATOM":
                self._natoms += 1
            else:
                break
        pdbfo.close()
        

    def setpdblist(self, pdblist):
        """Sets the list of PDBs.

        If a list is given, it must be a list of strings which are the PDBs
        to read in.  If a string is given, the string is expanded according
        to shell-globbing rules (*, ?, and [] expanded), the list is sorted,
        and the resulting list is used.  Note that if you don't have the PDBs
        in sort-order (eg. you have ... mol8.pdb, mol9.pdb, mol10.pdb),
        you've got to make a list yourself, put it in the right order, and
        pass the list to this function.
        """

        # first off, set self.pdblist to be the list of pdbs we get our
        # data from
        if isinstance(pdblist, str):
            pdblist = glob.glob(pdblist)
            pdblist.sort()
            self._pdblist = pdblist
        else:
            #This doesn't check to be sure that pdblist is a valid type.
            self._pdblist = pdblist
        self._initted = False
        self.init()
        
    def nextframe(self):
        """
        """
        self._framen += 1
        frame = numpy.zeros(shape=(self._natoms,3), dtype=numpy.float32)

        ### read through the PDB, 
        pdbfo = file(self._pdblist[self._framen], "r")
        while True:
            line = pdbfo.readline()
            if line[0:6] == "CRYST1":
                self._boxsize = numpy.array([float(line.split()[1]),
                                             float(line.split()[2]),
                                             float(line.split()[3])])
            if line[0:4] == "ATOM":   # we find the first atom line here, but
                                      #don't process it until the next loop.
                break
        # We will always keep atomn aligned with what line we have last read.
        atomn = 0
        try:
            while True:
                frame[atomn,0] = float(line[30:38])
                frame[atomn,1] = float(line[38:46])
                frame[atomn,2] = float(line[46:54])
                line = pdbfo.readline()
                atomn += 1
                if atomn+1 > self._natoms:
                    break    # let's be sure we get this right. when
                             # atomn+1 == natoms, we are on the last atom.  But
                             # we still have to go around one more time to be
                             # sure to process the current line, thus we use the
                             # inequality here.
        except:
            print "Exception at atom %s (starting at zero), line is: \n%s"%(atomn, line)
            raise

        self._frame = frame
        return self._frame


    def zero_frame(self):
        """Returns to frame zero.
        """
        self._framen = -1
        try:
            del self._frame
        except:
            pass
        
    
    def read_n_frames(self, number):
        """Reads `number' more frames, only printing out the last one.

        """
        self._framen += number-1
        return self.nextframe()

    #def __getattr__(self, attrname):
    #    """Wrapper for getting cord attributes.
    #
    #    This is a wrapper for getting things like self.nframes when you
    #    don't have to worry about setting them yourself.  Going on the
    #    hypothesis that most of the time these aren't changed, for
    #    anything that you haven't defined yourself, it will pass it
    #    through to self.cord.
    #
    #    Note that this could be bad in some cases!  But I'll take care
    #    of them when I find them.
    #    """
    #    return self.cord.__dict__[attrname]

