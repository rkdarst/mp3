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


class CordTinkerArc(mp3.cord.Cord):
    """Manages coordinates from a Tinker .arc files
    """
    initted = False
    def init(self):
        """Initilizes data about the PDBs.

        You don't need to use this, use self.setpdblist and this will
        be called automatically.

        The file can not have any other miscellaneous junk in it. 
        """
        # check to see if we have already been initilized, if so, don't
        # do it again.
        if self._initted == True:
            return None
        else:
            self._initted = True

        # number of atoms: it's the only thing on the first line
        self._natoms = int( file(self._arcfile).readline().strip() )

        # find number of frames and number of atoms
        # This is more tricky since it isn't encoded, we have to get it
        # from the length of the file.
        natoms = self._natoms
        number_of_lines = 0
        for line in file(self._arcfile, "r"):
            number_of_lines += 1
        lines_per_frame = natoms + 1
        #print "nlines", number_of_lines
        #print "lines_per_frame", lines_per_frame
        nframes = number_of_lines / lines_per_frame
        #print "nframes", nframes
        self._nframes = nframes
        self._framen = -1

        # we'll use a persistent file object for reading.
        # This may cause problems with pickling
        self._arc_fileobject = file(self._arcfile, "r")

    def __init__(self, arcfile):
        """Initilize the Tinker XYZ object.

        The argument is passed straight to setarc().
        """
        self.setarc(arcfile)

    def setarc(self, arcfile):
        """Sets the archive file to read.
        """

        # first off, set self.xyzlist to be the list of xyzs we get our
        # data from
        self._arcfile = arcfile
        self._initted = False
        self.init()
        
    def nextframe(self):
        """
        """
        self._framen += 1

        frame = numarray.zeros(shape=(self._natoms,3), type=numarray.Float32)

        ### read through the PDB, 
        arc_fo = self._arc_fileobject
        arc_fo.readline()   # bypass the number of atoms line
        for atomn in xrange(self._natoms):
            line = arc_fo.readline()
            frame[atomn,0] = float(line[12:23])
            frame[atomn,1] = float(line[24:35])
            frame[atomn,2] = float(line[36:47])
        self._frame = frame
        return self._frame


    def zero_frame(self):
        """Returns to frame zero.
        """
        self._framen = -1
        self._arc_fileobject.seek(0)
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
        return self.cord.__dict__[attrname]

