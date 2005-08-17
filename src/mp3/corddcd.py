import sys, string, struct
import mp3.log as log
log.debug('Loading corddcd.py')
import numarray
import mp3.cord
"""
This module contains definitions that pertain only to
single dcd files. Inherets from cord, which has
more general commands (like dcd printing)
"""



class CordDCD(mp3.cord.Cord):

    def __init__(self, dcd=None):
        """Create the CordDCD object.

        The keyword argument 'dcd' should be the filename to read.
        """
        # This is the initilization method for a dcd-cord object.  This
        # basically lets you set the dcd when you initilize the class.
        #
        if dcd != None:
            self.setdcd(dcd)
            self.init()

    def init(self):
        """Initilize by reading in the header information.

        You shouldn't need to directly call this anymore, unless you
        are doing something wierd.  It won't hurt to call this more
        than once.  It is now called automatically by setinput().
        """
        log.debug('--in corddcd.py, init()')
        if hasattr(self, '_initted') and self._initted == True:
            return None
        self._get_header()
        self._initted = True  #make a record saying that we've been set up already
        self._framen = -1   #THE FIRST FRAME IS FRAME 0--- just like python
                           # BUT!!!!!!!!
                           # the actual DCD does not contain the initial
                           # coordinates.  So maybe this should be "one" here.
                           # That will change lots of stuff.  I'll leave it
                           # like this for now.  

        self._frame = numarray.zeros(shape=(self._natoms,3), type=numarray.Float32)


    ##def load(self):
    ##    """load the entire dcd into memory"""
    ##    thelog.debug('--In dcd.py, load()')
    ##    #self.init()   #if it's already initted, it'll catch itself. 
    ##    self.load_vectors()

    def nextframe(self):
        """Loads the next frame into self.
        loads the next frame in self.frame . Hides details about wheter
        or not the dcd has been loaded into memory
        """
        #log.debug('--In dcd.py, nextframe()')
        #if self.initted == 0:
        #    thelog.critical("init it first! I won't do it here to save time...")
        #    return
        #thelog.info('make a better solution for this...')
        if self._framen >= self._nframes-1:  #if there are 10 frames total, framen can not be
                                           #more than 9. So stoup if framen >= nframes -1
             log.critical("tried to read more frames than there were. BAILING OUT!!")
             return None
        ##if hasattr(self, 'data'):
        ##    self.framen += 1
        ##    self.frame = self.data[self.framen]
        ##    return self.frame
        self._get_next_frame()
        self._framen += 1
        return self._frame

    def read_n_frames(self, ntoread):
        """Advance the current frame by the given number.

        
        Reads n frames from the dcd, only returning ( and saving
        in self.frame) the last one.  If called with an argument of one,
        it behaves exactly like nextframe()

        Note that it might just skip the data for speed issues, it need not
        concern you.
        """
        #log.debug('--In dcd.py, read_n_frames()')
        # note that this has to accomodate both the cases where you
        # have already loaded it all into memory and where you are
        # reading it frame-by-frame. You could either call self.nextframe()
        # the right number of times, or look at which case it is and optimize
        # for the skipping. I'll do the optimized method (though I'm sure
        # that it makes very little difference. _very_. )
        # 
        # I'm also violating a fundamental rule of python by not combining
        # this function with the one above, but I'll ignore that for now.
        # 
        # Update: if we did it by calling self.next_frame() repeatedley,
        # we could make this part of the generalized dcd class
        ntoread = int(ntoread)
        if self._framen + ntoread > self._nframes - 1:   #can't read more than we have
            log.critical("tried to read more frames than there are!")

#I'll advance the dcd file by the right amount
        ntoskip = ntoread-1
        self._fo.seek( (self._natoms*4 + 8 )*3*ntoskip, 1 )  #seek, whence=1 ==> seek from current position
        self._framen += ntoskip
        self._get_next_frame()
        self._framen += 1
        return self._frame

    def zero_frame(self):
        """makes it where the next frame returned is 0. This function itself
        does not return anything.
        """
        log.debug('---In dcd.py, zero_frame(). self.framen=%s'%self.framen())
        ##if self.loadedinmem == 1:
        ##    thelog.debug("dcd.py, zero_frame(), it's loaded in mem so just resetting framen")
        ##    self.framen = -1    #we don't have to deal with any file rewinding here
        ##else:
        log.debug("dcd.py, zero_frame(), resetting file object and framen")
        self._framen = -1
            #self.fo.seek(0)     #we have to seek past the header somehow.
            #(combined below)    #but actual header size varries with the length
                                 #of the title. So I'll do a quick calculation to see
                                 #how much we have to skip. It's ugly, but it works.
                                 #
                                 #84 + 8 = 92 bytes for the first header
                                 #len(title) + 12 for the lenght of the title header
                                 #12 bytes for the lenght of the atom num header
                                 #total = 116 + len(title)
        self._fo.seek(116 + len(self._title) ) #seek from beginning of the file

    def setdcd(self, filename):
        """Set the filename of this dcd, from a string

        sets the filename (string, not file object) that we will be
        reading from.  Now, it also .init()s the dcd object.
        """
        log.debug('---in dcd.py, setfo()')
        self._fo = file(filename,"r")
        self.init()
        

##
##  Reading-type commands.
##
##
    def _get_header(self):
        """Reads in the header of the dcd.

        Not for public use.
        Wrapper method. Calls methods to read in all the headers.
        """
        log.debug('---In dcd.py, get_header()')
        self._parse_header()
        self._parse_title()
        self._parse_atoms()
    
        log.debug("Parsing: %s"% self._title)    #print out some useful information
        for i in range(0,len(self._title),80):
            log.debug(self._title[i:i+80])

        if self._nframes*self._dcdfreq != self._ntsteps:
            log.warn("error-- the wierd ntsteps frame is not what I think it should be!")
    
    ####end main funcition####################
    def _parse_header(self):
        """Parses in the first header.

        Not for public use
        """
        log.debug('---In dcd.py, parse_header()')
        #process the first header block

        header1 = self._fo.read(92)
        header1_format=\
        "i---cccci---i---i---i---xxxxxxxxxxxxxxxxxxxxf---i---i---xxxxxxxxxxxxxxxxxxxxxxxxxxxxi---i---"
        # |1  |5   |10  |15  |20  |25  |30  |35  |40  |45  |50  |55  |60  |65  |70  |75  |80  |85  |90
        #|header size=84     |nframes*tstep          |tstep_size                             |charm_ver
        #    |CORD=has coordinates                       |block_a                                |header_size=84
        #        |nframes                                    |block_b
        #            |starting timestep
        #                |timestep between coord sets                                                  
        header1_format = string.replace(header1_format, "-", "")
        header1 = struct.unpack(header1_format, header1)
        header1_size1, c1, c2, c3, c4, self._nframes, self._firsttstep, self._dcdfreq, self._ntsteps, self._tstep_size, self._block_a, self._block_b, self._charm_v, header1_size2 = header1  #unpack the tuple header1
    
    
        self._dcdtype = "".join((c1,c2,c3,c4))   #get the data-type field. I it should always be cord...
        if header1_size1 != 84 or header1_size2 !=84:
            log.error("error-- header size fields not correct (should be 84)\n")
        if self._block_a != 0 or self._block_b != 0:
            log.error("I've found some charmm wierdness-- the alignment could be all wrong-- no corrections will be made!")
    
    #########################
    def _parse_title(self):
        """Parses in the title.

        Not for public use.
        """
        log.debug('---In dcd.py, parse_title()')
        ##
        ##Parse the title field
        ##
        #title field format: int that specifies the block size = 80n + 4
        #                    int that specifies the number of 80-char title blocks = x
        #                    x * 80 char title blocks
        #                      ...
        #                    int that specifies the block size = 80n +4
        title = ""
        blocksize, n = struct.unpack("ii", self._fo.read(8))   #read in 8 bytes, convert to ints, store in blocksize, n
        if blocksize != 80*n+4:
            log.error("beginningblock size in title block is wrong\n")
        for i in range(0,n):
            title = title + self._fo.read(80)          # read the next title string onto the current one.
        if (blocksize, ) != struct.unpack("i", self._fo.read(4)):
            log.error("ending blocksize on the title block is not correct !\n")
        self._title = title
    
    ######################################
    
    ###########################
    def _parse_atoms(self):
        """Parses in the natoms record.

        Not for public use.
        """
        ##
        ##Parse the number of atoms
        ##
        #format: 3 ints, first one is blocksize (must equal 4), second is natoms, third is blocksize (must equal 4)
        log.debug("---in dcd.py, parse_atoms()")
        blocksize, self._natoms, blocksize2 = struct.unpack("iii", self._fo.read(12))
        if blocksize != 4 or blocksize2 != 4:
            log.error("blocksizes in the number of atoms record is broken\n")
        
    
    ##########################
    
    ##def load_vectors(self):
    ##    thelog.debug('--In dcd.py, load_vectors()')
    ##    ##
    ##    ##Parse the coordinate/velocity fields
    ##    ##
    ##    #I won't go into the format again.
    ##
    ##    #  here is the plan. Make a array that is big enough for all of the data that we will get,
    ##    # then assign the right values to slices of the array
    ##
    ##    self.data = numarray.zeros(shape=(self.nframes,self.natoms,3), type=numarray.Float32)     #define the data array
    ##
    ##    for framen in range(0,self.nframes):
    ##        #for coordinate in (0,1,2):
    ##        self.fo.read(4)             #to save time, we won't error check this ! (or should we ?
    ##        self.data[framen,:,0] = numarray.fromfile(self.fo, numarray.Float32, shape=(self.natoms,))
    ##        self.fo.read(8)
    ##        self.data[framen,:,1] = numarray.fromfile(self.fo, numarray.Float32, shape=(self.natoms,))
    ##        self.fo.read(8)
    ##        self.data[framen,:,2] = numarray.fromfile(self.fo, numarray.Float32, shape=(self.natoms,))
    ##        self.fo.read(4)
           
    ###########################
    
    def _get_next_frame(self):
        """Read in next frame, store it as self.frame.
        
        Not for public use.
        This reads in the next frame from the dcd's file object (self.fo), and
        stores it as self.frame.  No return value.  Does not advance self.framen.
        Look at nextframe() instead.
        """
        #log.debug('---In dcd.py, get_next_frame(), self.framen=%s'%self.framen())
        self._fo.read(4)             #to save time, we won't error check this ! (or should we ?
        self._frame[:,0]= numarray.fromfile(self._fo, numarray.Float32, shape=(self._natoms,))
        self._fo.read(8)
        self._frame[:,1] = numarray.fromfile(self._fo, numarray.Float32, shape=(self._natoms,))
        self._fo.read(8)
        self._frame[:,2] = numarray.fromfile(self._fo, numarray.Float32, shape=(self._natoms,))
        self._fo.read(4)
