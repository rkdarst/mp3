# This module serves as an interface to PDBs.  It allows you to
# access the coordinates stored in PDBs.
#
#

import logging ; thelog = logging.getLogger('mp3')
thelog.debug('Loading pdbcord.py')
import glob
import numarray
import numarray
import mp3.cord


class CordTinkerXYZ(mp3.cord.Cord):
    """Manages coordinates from a series of Tinker .xyz files
    """
    initted = False
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
        self._nframes = len(self._xyzlist)
        self._framen = -1
        # number of atoms: it's the only thing on the first line
        self._natoms = int( file(self._xyzlist[0]).readline().strip() )

    def __init__(self, xyzlist):
        """Initilize the Tinker XYZ object.

        The argument is passed straight to setxyzlist().
        """
        self.setxyzlist(xyzlist)

    def setxyzlist(self, xyzlist):
        """Sets the list of .xyz's.

        If a list is given, it must be a list of strings which are the PDBs
        to read in.  If a string is given, the string is expanded according
        to shell-globbing rules (*, ?, and [] expanded), the list is sorted,
        and the resulting list is used.  Note that if you don't have the PDBs
        in sort-order (eg. you have ... mol8.xyz, mol9.xyz, mol10.xyz),
        you've got to make a list yourself, put it in the right order, and
        pass the list to this function.
        """

        # first off, set self.xyzlist to be the list of xyzs we get our
        # data from
        if type(xyzlist) == str:
            xyzlist = glob.glob(xyzlist)
            xyzlist.sort()
            self._xyzlist = xyzlist
        else:
            #This doesn't check to be sure that xyzlist is a valid type.
            self._xyzlist = xyzlist
        self._initted = False
        self.init()
        
    def nextframe(self):
        """
        """
        self._framen += 1
        frame = numarray.zeros(shape=(self._natoms,3), type=numarray.Float32)

        ### read through the PDB, 
        xyzfo = file(self._xyzlist[self._framen], "r")
        xyzfo.readline()   # bypass the number of atoms line
        for atomn in xrange(self._natoms):
            line = xyzfo.readline()
            frame[atomn,0] = float(line[12:23])
            frame[atomn,1] = float(line[24:35])
            frame[atomn,2] = float(line[36:47])
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
        #print attrname
        #raise
        print attrname
        raise
        return self.cord.__dict__[attrname]

