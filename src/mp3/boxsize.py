# Various miscellaneous tools for dealing with varying box sizes. 
#
# A box size is a tuple with elements 0,1,2 representing the cartesian
# x,y,z axes lenghts.  (Well, it doesn't have to be a tuple but
# boxsize[i] for i in [1,2,3] should return the right things.
# 
# Box Size protocol:
#   __call__(), returns the new box size (advances frame)
#   framen(),   which frame number it is on, starting from zero.
#   boxsize     attribute which contains the boxsize tuple

class BoxStatic:
    """A box size that never changes.

    
    """
    def __init__(self, boxsize):
        """Initilize and set the static box size.

        The keyword argument 'boxsize' is a tuple which is the box
        size.
        """
        self._boxsize = size
        self.framen = -1
    def __call__(self):
        self.framen += 1
        return self._boxsize

class BoxXST:
    """Interface to xst (extended system trajectory) files

    This class allows you to get a handle on xst files, and extract
    the frame sizes from them. Use them thus:

    box = xstbox("filename.xst")   #instantiate the object
      # this does:
      # -stores the filename
      # -sets sets box.boxsize to be the size of the zeroth frame
      #                      NOT the first frame of a DCD...
      #                      xst's have the initial coordinates, dcd's don't
      # - sets box.framen = 0  (since this is the zeroth frame...)

    box()   # returns the next frame's size (x,y,z).
            #also stores it in box.boxsize and advances box.framen
               # 'box' is a class that is callable 

    USEFUL METHODS/ATTRIBUTES:
    self.framen  --  the framenumber.  zero=initial cords, not in DCD
    self.boxsize  -- last called frame

    
    """
    def __init__(self, xst=None):
        """Sets it up.

        This method takes an optional argument, which is the xst filename.
        Giving the optional argument replaces the need to call .setxst()
        """
        if xst != None:
            self.setxst(xst)

    def setxst(self, xstfile):
        """Sets the underlying xstfile.

        Takes a string, and sets up what where the xst data comes from.

        
        """
        xst = file(xstfile,"r")
        xst.readline()
        xst.readline()
        self.xst = xst
        self.framen = -1
        self()
    def __call__(self):
        line = self.xst.readline().split()
        self.boxsize = (float(line[1]), float(line[5]), float(line[9]))
        self.framen += 1
        return self.boxsize
Xstbox = BoxXST
    
class BoxMultiXST:
    """Return box sizes when you have multiple XST files.
    """
    def __init__(self, xstlist):
        """

        The keyword argument 'xstlist' is the list of XST files to read from.
        """
        print "MULTI-XST BOXSE ARE EXPERIMENTAL"
        xstobjects = []
        lengths    = []
        for xst in xstlist:
            length = len(file(xst),"r").readlines() -2
            #lenght is the number of lines of data in the XST
            #lenght-1 is the actual number of lines in it.
            lengths.append(length-1)
            xstobjects.append(BoxXST(xst=xst))
        self._xstobjects = hi
        self._lenghts = lengths
        self._index = 0
    def __call__(self):
        i = self._index
        if self._xstobjects[i].framen >= self._lenghts[i]:
            self._index += 1
            return self()
        else:
            return self.xstobjects[i]()
        
    
