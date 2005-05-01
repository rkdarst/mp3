

class Xstbox:
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
    

        
