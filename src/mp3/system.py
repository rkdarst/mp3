import logging
thelog = logging.getLogger('mp3')
import labels #as mp3labels
import numpy

class System(labels.Labels):
    """Object which combines coordinates and atom labels.

    This object stores coordinates at self.cord, which are regular
    coordinate objects.  Label data is stored at self.data.

    Label data is accessed by self.data['labelname'][atomnum].  It is
    generally more efficient to do data['labelname'][atomnum] than
    data[atomnum][labelname].  Label data is implemented by the class
    mp3.labels.Labels.  The System class is a subclass of Labels, so
    all Labels methods are called directly on the system.  For
    example, you would do System.getfrompdb().  See `pydoc
    mp3.labels.Labels` for more information on using labels.
    """

    def __init__(self, psf=None, cord=None, pdb=None):
        """
        """
        self.labels = self
        if psf != None:
            self.getfrompsf(psf)
        if pdb != None:
            self.getfrompdb(pdb)
        if cord != None:
            self.cord = cord

    #def __getattr__(self, name):
    #    """Magic method to ease working with systems
    #
    #    This allows the following properties to be accessed from system objects.
    #    Otherwise, the user must use system.labels, or system.cord.  Note that
    #    if system.cord does not exist, exceptions will be raised.
    #    """
    #    # So the point of this is to let you easily access some key
    #    # things in self.labels and self.cord easily
    #
    #    cordattr = ("frame","nextframe","framen", "natoms")
    #    labelattr = ("findatoms","lab", "findrange")
    #    try:
    #        if name in cordattr:
    #            return getattr(self.cord, name)
    #        if name in labelattr:
    #            if name == "lab":
    #                return self.labels.data
    #            return getattr(self.labels, name)
    #    except:
    #        raise AttributeError
    #    raise AttributeError
    def natoms(self):
        """Number of atoms.

        First attempts to get number of atoms from labels which have
        been read, if that is not avaliable then it will attempt to
        get it from the number of atoms in a coordinate file.
        """
        if hasattr(self, "_natoms"):
            return self._natoms
        if hasattr(self, "cord"):
            return self.cord.natoms()
    #def __set_natoms(self, natoms): self.natoms = natoms
    #natoms = property(__get_natoms, __set_natoms, None, """Get number of atoms in the system.
    #
    #This searches self for the attribute _natoms (this is set when
    #PSFs or other data are read), if that is not found, it will then
    #search self.cord for natoms.
    #""")

    def setcord(self, cord):
        """Set the location to get input cordinates.

        The input coordinates are found at self.cord
        """
        self.cord = cord

#
#  Writing PDBs
#

    def writepdbseries(self, prefix, atomlist=None):
        """Writes out an entire pdb sequence.

        Pass this method the prefix of your PDBs.  It will write out
        all frames in the coordinate set in self.cord, using field
        information from self.data .

        Prefix should be a string. Output is "prefixXXXXX.pdb", where
        XXXXX is the number of the frame, starting from zero.

        You can use the `atomlist` option to write out only a subset
        of the atoms.

        This method will call cord.nextframe() before it starts printing,
        so do not call self.cord.nextframe() yourself-- this will break
        things (it will try to read one too many frames).
        """
        thelog.debug("system.py, printpdbseq()")
        for framen in range(0,self.cord.nframes()):
            filename = prefix + ("%.4d" % framen) + ".pdb"

            self.cord.nextframe()
            self._pdbframe(filename, atomlist=atomlist)

    def writepdb(self, name, atomlist=None):
        """Writes a single pdb, from the current frame.

        You can use the `atomlist` option to write out only a subset
        of the atoms.

        Writes the current frame (self.cord.frame) to the given
        filename (not just a prefix).  nextframe() is not called.       
        """
        self._pdbframe(name, atomlist=atomlist)

# 
# Stuff for printing PBDs
#     
    def _pdbframe(self, name, atomlist=None):
        """Writes the current frame to a PDB file.

        Use <system>.writepdb() instead.
        """
        self.output = file(name,"a")
        self.output.write("MODEL\n")        
        if not hasattr(self, "atomlist"):  #self.use_atom_subset == False:
            for atomn in range(0, self.cord.natoms()):
                self._pdbline(atomn, atomn)
        else: # we are not using a subset  # self.use_atom_subset == True:
            counter = 0  # it is increased by one in the pdbline function
            for atomn in self.atomlist:
                self._pdbline(atomn, counter)
                counter += 1
        self.output.write("ENDMDL\n")
        self.output.close()
        del self.output

     
    def _pdbline(self, atomn, fakeatomn):
        """Writes a line from a pdb.
        
        non-public method. Writes the given pdb line to the file object
        attached.

        This method now requires two arguments. The first one is the atom in
        the DCD/labels that we want to print, the second one is what number we
        want to print it as. The reason for this is to make it easier to
        print out a subset of atoms-- we can treat everything equally.

        """
    
        #ATOM      1  N   MET X   1      27.340  24.430   2.614  1.00  9.67      1UBQ
        #|0   |5   |10  |15  |20  |25  |30  |35  |40  |45  |50  |55  |60  |65  |70  |75  |80
        #|1  |5   |10  |15  |20  |25  |30  |35  |40  |45  |50  |55  |60  |65  |70  |75  |80
        #s-----i---- s--- s-- si---c   f-------f-------f-------f-----f-----      s---s-s-
        #ATOM     23  CG1 ILE X   3      27.810  28.748   8.999  1.00  4.72      1UBQ
        #ATOM     44  HA   SER     3     -17.929  -9.853  12.085  0.00  0.00      APO
        #ATOM                 |chainid         |y                                |segname
        #      |atomn          |resnum                 |z                            |element
        #            |atomtype     |icode?                     |occupancy              |charge
        #                 |resname     |x                            |tempfactor
       #"s-----

        output = self.output  #=self.output

        ##  output.write( "ATOM  " )
        ##  output.write( "%5d " % (fakeatomn +1) )  #atomn
        ##  output.write( "%-4.4s " % self.data['atomname'][atomn] ) 
        ##  # (above)
        ##  # it seems that a numarray.recarray strips any tailing whitespace.
        ##  # This natually will break the spacing of the output if you try to
        ##  # write it right-justified.  So we can't do that.  So I will make it
        ##  # left-justified until something better comes up.
        ##  # XXX fixed in numpy-- I should come back and reanalyze this write-out.
        ##  
        ##  output.write( "%-4.4s" % self.data['resname'][atomn] )
        ##  #this is the field which should be three chars, but I illegally increase it to four
        ##  output.write(" "  ) #chainid 
        ##  output.write( "%4d" % self.data['resnum'][atomn] )
        ##  output.write(" "  )#icode
        ##  output.write("   ")
        ##  output.write( "%8.3f" % self.cord.frame()[atomn,0] )    #use "%8.3f" for proper pdb format  #xcord
        ##  output.write( "%8.3f" % self.cord.frame()[atomn,1] ) #ycord
        ##  output.write( "%8.3f" % self.cord.frame()[atomn,2] ) #zcord
        ##  output.write( "%6.2f" % self.data['occupancy'][atomn] )     #occupancy
        ##  output.write( "%6.2f      " % self.data['tempfactor'][atomn] )    #tempfactor
        ##  output.write( "%-4.4s" % self.data['segname'][atomn] ) #segname
        ##  output.write( "%2.2s"% self.data['element'][atomn] )  #element
        ##  
        ##  output.write( "  \n" )   #charge

        data = self.data
        atomname = "%-4.4s " % data['atomname'][atomn]
        #atomname = data['atomid'][atomn]
        #if len(atomname) == 4:
        #    atomname = "%4.4s " % atomname
        #else:
        #    atomname = " %-3.3s " % atomname
        segname = "%-4.4s" % self.data['segname'][atomn]
        #segname = data['segname'][atomn]
        #if len(segname) == 4:
        #    segname = "%-4.4s" % segname
        #else:
        #    segname = " %3.3s" % self.data['segname'][atomn]
        line = ["ATOM  ",
                "%5d "    % (fakeatomn +1),              #atomn
                atomname,                                #atomname  
                "%-4.4s"  % data['resname'][atomn],      #resname
                " ",                                     #chainid
                "%4d"     % data['resnum'][atomn],       #resnum
                " ",                                     #icode
                "   ",
                "%8.3f"   % self.cord.frame()[atomn,0],  #xstr  #use "%8.3f" for proper pdb format
                "%8.3f"   % self.cord.frame()[atomn,1],  #ystr
                "%8.3f"   % self.cord.frame()[atomn,2],  #zstr
                "%6.2f"   % data['occupancy'][atomn],    #occupancy
                "%6.2f      " % data['tempfactor'][atomn], #tempfactor
                segname,                                 # segname
                "%2.2s"   % data['element'][atomn],      #element
                "  \n",    #charge
                ]
        output.write("".join(line))

    def _pdbline_broke1(self, atomn, fakeatomn):
        """Just like pdbline.  non-public.

        So what this does differently is it aligns some things slightly
        differently.  The reason it exists is to accomodate a broken program.
        This: http://www.rcsb.org/pdb/docs/format/pdbguide2.2/part_13.html
        says that alignment should not matter.

        To use it, do this:
        mysystem = mp3.system()
        mysystem._pdbline = mysystem._pdbline_broke1
        """
    
        #ATOM      1  N   MET X   1      27.340  24.430   2.614  1.00  9.67      1UBQ
        #|0   |5   |10  |15  |20  |25  |30  |35  |40  |45  |50  |55  |60  |65  |70  |75  |80
        #|1  |5   |10  |15  |20  |25  |30  |35  |40  |45  |50  |55  |60  |65  |70  |75  |80
        #s-----i---- s--- s-- si---c   f-------f-------f-------f-----f-----      s---s-s-
        #ATOM     23  CG1 ILE X   3      27.810  28.748   8.999  1.00  4.72      1UBQ
        #ATOM     44  HA   SER     3     -17.929  -9.853  12.085  0.00  0.00      APO
        #ATOM                 |chainid         |y                                |segname
        #      |atomn          |resnum                 |z                            |element
        #            |atomtype     |icode?                     |occupancy              |charge
        #                 |resname     |x                            |tempfactor
       #"s-----

        output = self.output  #=self.output

        data = self.data
        atomname = data['atomid'][atomn]
        if len(atomname) == 4:
            atomname = "%4.4s " % atomname
        else:
            atomname = " %-3.3s " % atomname
        segname = data['segname'][atomn]
        if len(segname) == 4:
            segname = "%-4.4s" % segname
        else:
            segname = " %3.3s" % self.data['segname'][atomn]
        line = ["ATOM  ",
                "%5d "    % (fakeatomn +1),              #atomn
                atomname,                                #atomname  
                "%-4.4s"  % data['resname'][atomn],      #resname
                " ",                                     #chainid
                "%4d"     % data['resnum'][atomn],       #resnum
                " ",                                     #icode
                "   ",
                "%8.3f"   % self.cord.frame()[atomn,0],  #xstr  #use "%8.3f" for proper pdb format
                "%8.3f"   % self.cord.frame()[atomn,1],  #ystr
                "%8.3f"   % self.cord.frame()[atomn,2],  #zstr
                "%6.2f"   % data['occupancy'][atomn],    #occupancy
                "%6.2f      " % data['tempfactor'][atomn], #tempfactor
                segname,                                 # segname
                "%2.2s"   % data['element'][atomn],      #element
                "  \n",    #charge
                ]
        output.write("".join(line))

    def writetinkerxyz(self, name):

        cords = self.cord.frame()
            
        fo = file(name, "w")
        data = self.data
        fo.write("%d\n"%self.natoms())

#   11  HC    -0.023194    1.941312    0.866941     6     9
      
        for i in range(self.natoms()):
            fo.write("%6d  %-4.4s%11.6f %11.6f %11.6f %5.5d%s"%(i+1,
                                                                data["atomname"][i],
                                                                cords[i,0],
                                                                cords[i,1],
                                                                cords[i,02],
                                                                data["atomtypenum"][i],
                                                                self.labels._tinkerbondlist[i]  ))

    def center_of_mass(self, atomlist=None, weights=None):
        """Calculate the center of mass of the system.

        If 'atomlist' is None, use all atoms to calculate center of
        mass.  Otherwise, atomlist must be a list containing the atoms
        to use in calculating the center of mass.

        If 'weighting' is not None, do not weight by masses, but
        instead weight by the weighting array, which should be an
        array with len(atomlist) atoms (Note: it is not sliced, if
        atomlist is passed, weighting must already have the proper
        shape)
        """
        # get atomlist if not passed, defaulting to all atoms.
        if atomlist is None:
            atomlist = range(self.natoms())

        # get weighting if not passed, defaulting to masses from the
        # labels.
        if weights is None:
            weights = self.data["mass"][atomlist]
        # get positions for our atoms
        #   I tested this, it won't alter the original frame.
        frame = self.cord.frame()[atomlist]
        frame.transpose()
        # COM = sum(r*m) / sum(m)
        COM = numpy.sum(frame*weights, axis=1) / sum(weights)

        # This code has been tested against VMD.  This is what we got.
        # It isn't exactly on, but close enough to blame it on
        # numerical error.  
        #
        # # file is 'files/apombC_278Kc_1-2ab_WRAP_last1ns.dcd'
        # # all atoms
        # vmd > measure center $prot weight mass
        # 0.058566596359 0.0782859697938 -0.0550067648292
        # >>> COM
        # array([-0.07915176,  0.07452537, -0.04425772], type=Float32)
        #
        # # just a protein
        # vmd > measure center $prot weight mass
        # 0.782196938992 12.0635070801 -6.97437620163
        # >>> COM
        # array([  0.78218687,  12.06335163,  -6.97428656], type=Float32)
        return COM
