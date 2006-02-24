import logging
thelog = logging.getLogger('mp3')
import labels #as mp3labels

import numarray

class System:

    def __init__(self, psf=None, cord=None, pdb=None):
        """
        """
        self.labels = labels.Labels()
        if psf != None:
            self.labels.getfrompsf(psf)
        if pdb != None:
            self.labels.getfrompdb(pdb)
        if cord != None:
            self.cord = cord
            

    def __getattr__(self, name):
        """Magic method to ease working with systems
    
        This allows the following properties to be accessed from system objects.
        Otherwise, the user must use system.labels, or system.cord.  Note that
        if system.cord does not exist, exceptions will be raised.
        """
        # So the point of this is to let you easily access some key
        # things in self.labels and self.cord easily

        cordattr = ("frame","nextframe","framen", "natoms")
        labelattr = ("findatoms","lab", "findrange")
        try:
            if name in cordattr:
                return getattr(self.cord, name)
            if name in labelattr:
                if name == "lab":
                    return self.labels.data
                return getattr(self.labels, name)
        except:
            raise AttributeError
        raise AttributeError

    def setcord(self, cord):
        """Set the location to get input cordinates.
        """
        self.cord = cord

#
#  Writing PDBs
#

    def writepdbseries(self,prefix):
        """Writes out an entire pdb sequence.

        Pass this method the prefix of your PDBs.  It will write out all
        frames in the coordinate set in self.cord, using field
        information from self.labels .

        Prefix should be a string. Output is "prefixXXXXX.pdb", where
        XXXXX is the number of the frame, starting from zero. You can use
        .atoms_to_use() to print out only a subset of the atoms.

        This method will call cord.nextframe() before it starts printing,
        so do not call self.cord.nextframe() yourself-- this will break
        things (it will try to read one too many frames).
        """
        self._writepdbseq(prefix)

    def writepdb(self, name):
        """Writes a single pdb, from the current frame.

        Writes the current frame (self.cord.frame) to the given
        filename (not just a prefix).  If a atomlist has been defined
        using .set_atoms_to_print(), it will only print out lines
        corresponding to those atoms        
        """
        self._pdbframe(name)

# 
# Stuff for printing PBDs
# 
    def atoms_to_use(self, atomlist=None):
        """Set up for printing only a subset of all atoms in the PDBs

        If this method has been called and given an atomlist, the list
        will be stored and any PDBs printed will only include those atoms.

        If this method is called without an argument or with None as
        the argument, it will reset the state so that all atoms will
        be included.        
        """
        if atomlist != None:
            self.atomlist = atomlist
        else:
            if hasattr(self, "atomlist"):
                del self.atomlist

    def _writepdbseq(self, prefix):
        """Writes out an entire pdb sequence. (non-public)

        Use .writepdbseries() instead.  That's a wrapper that
        won't change, this function might
    
        prints all pdb frames in a given system, with the given prefix
        (prefix is a string).

        This method will call self.cord.nextframe() before it writes out it's
        first frame!

        Not a public method. 
        """
        thelog.debug("system.py, printpdbseq()")
        for framen in range(0,self.cord.nframes()):
            filename = prefix + ("%.4d" % framen) + ".pdb"

            self.cord.nextframe()
            self._pdbframe(filename)
    
    
    def _pdbframe(self, name):
        """Writes the current frame to a PDB file.

        Use <system>.writepdb() instead.
        """
        self.output = file(name,"w")
        if not hasattr(self, "atomlist"):  #self.use_atom_subset == False:
            for atomn in range(0, self.cord.natoms()):
                self._pdbline(atomn, atomn)
        else: # we are not using a subset  # self.use_atom_subset == True:
            counter = 0  # it is increased by one in the pdbline function
            for atomn in self.atomlist:
                self._pdbline(atomn, counter)
                counter += 1
        self.output.write("END\n")
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

        output.write( "ATOM  " )
        output.write( "%5d " % (fakeatomn +1) )  #atomn

        output.write( "%-4.4s " % self.labels.data.field('atomname')[atomn] ) 
        # (above)
        # it seems that a numarray.recarray strips any tailing whitespace.
        # This natually will break the spacing of the output if you try to
        # write it right-justified.  So we can't do that.  So I will make it
        # left-justified until something better comes up.

        output.write( "%-4.4s" % self.labels.data.field('resname')[atomn] )
        #this is the field which should be three chars, but I illegally increase it to four
        
        output.write(" "  ) #chainid 
        output.write( "%4d" % self.labels.data.field('resnum')[atomn] )
        output.write(" "  )#icode
        output.write("   ")
        output.write( "%8.3f" % self.cord.frame()[atomn,0] )    #use "%8.3f" for proper pdb format  #xcord
        output.write( "%8.3f" % self.cord.frame()[atomn,1] ) #ycord
        output.write( "%8.3f" % self.cord.frame()[atomn,2] ) #zcord
        output.write( "%6.2f" % self.labels.data.field('occupancy')[atomn] )     #occupancy
        output.write( "%6.2f      " % self.labels.data.field('tempfactor')[atomn] )    #tempfactor
        output.write( "%-4.4s" % self.labels.data.field('segname')[atomn] ) #segname
        output.write( "%2.2s"% self.labels.data.field('element')[atomn] )  #element
        output.write( "  \n" )   #charge

#output.write( ( "ATOM  %s %s %s %s%s%s   %s%s%s%s%s      %s%s%s\n" % (atomserial, atomtype,resname,chainid,resnum,icode,xstr,ystr,zstr,occupancy,tempfactor,segname,element,charge))


#        atomserial = "%5d" % (atomn +1)
#        atomtype = "%4.4s" % self.labels.data.field('atomid')[atomn] #used to be atomtype
#        resname = "%3.3s" % self.labels.data.field('resname')[atomn]    #take the slice to ensure that it is 4 characters long
#        chainid = " "
#        resnum = "%4d" % self.labels.data.field('resnum')[atomn] 
#        icode = " "
#        xstr = "%8.4f" % self.cord.frame[atomn,0]    #use "%8.3f" for proper pdb format
#        ystr = "%8.4f" % self.cord.frame[atomn,1]
#        zstr = "%8.4f" % self.cord.frame[atomn,2]
#        occupancy = "%6.2f" % 0
#        tempfactor = "%6.2f" % 0
#        segname = "%4.4s" % self.labels.data.field('segname')[atomn]
#        element = "  "
#        charge = "  "
#        output = ( "ATOM  %s %s %s %s%s%s   %s%s%s%s%s      %s%s%s\n" % (atomserial, atomtype,resname,chainid,resnum,icode,xstr,ystr,zstr,occupancy,tempfactor,segname,element,charge))
#
        #return output 
        #self.output.write(output)

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

        output.write( "ATOM  " )
        output.write( "%5d " % (fakeatomn +1) )  #atomn

        if len(self.labels.data.field('atomid')[atomn]) == 4:
            output.write( "%4.4s " % self.labels.data.field('atomname')[atomn] ) #used to be atomtype  #atomtype
        else:
            output.write( " %-3.3s " % self.labels.data.field('atomtype')[atomn] ) #used to be atomtype  #atomtype



        output.write( "%-4.4s" % self.labels.data.field('resname')[atomn] )   #take the slice to ensure that it is 4 characters long  #resname
        output.write(" "  ) #chainid 
        output.write( "%4d" % self.labels.data.field('resnum')[atomn] ) #resnum
        output.write(" "  )#icode
        output.write("   ")
        output.write( "%8.3f" % self.cord.frame()[atomn,0] )    #use "%8.3f" for proper pdb format  #xstr
        output.write( "%8.3f" % self.cord.frame()[atomn,1] ) #ystr
        output.write( "%8.3f" % self.cord.frame()[atomn,2] ) #zstr
        output.write( "%6.2f" % self.labels.data.field('occupancy')[atomn] )     #occupancy
        output.write( "%6.2f      " % self.labels.data.field('tempfactor')[atomn] )    #tempfactor

        if len(self.labels.data.field('segname')[atomn]) == 4:
            output.write( "%-4.4s" % self.labels.data.field('segname')[atomn] ) #segname
        else:
            output.write( " %3.3s" % self.labels.data.field('segname')[atomn] ) #segname


        output.write( "  " )  #element
        output.write( "  \n" )   #charge

#output.write( ( "ATOM  %s %s %s %s%s%s   %s%s%s%s%s      %s%s%s\n" % (atomserial, atomtype,resname,chainid,resnum,icode,xstr,ystr,zstr,occupancy,tempfactor,segname,element,charge))


#        atomserial = "%5d" % (atomn +1)
#        atomtype = "%4.4s" % self.labels.data.field('atomid')[atomn] #used to be atomtype
#        resname = "%3.3s" % self.labels.data.field('resname')[atomn]    #take the slice to ensure that it is 4 characters long
#        chainid = " "
#        resnum = "%4d" % self.labels.data.field('resnum')[atomn] 
#        icode = " "
#        xstr = "%8.4f" % self.cord.frame[atomn,0]    #use "%8.3f" for proper pdb format
#        ystr = "%8.4f" % self.cord.frame[atomn,1]
#        zstr = "%8.4f" % self.cord.frame[atomn,2]
#        occupancy = "%6.2f" % 0
#        tempfactor = "%6.2f" % 0
#        segname = "%4.4s" % self.labels.data.field('segname')[atomn]
#        element = "  "
#        charge = "  "
#        output = ( "ATOM  %s %s %s %s%s%s   %s%s%s%s%s      %s%s%s\n" % (atomserial, atomtype,resname,chainid,resnum,icode,xstr,ystr,zstr,occupancy,tempfactor,segname,element,charge))
#
        #return output 
        #self.output.write(output)


    def writetinkerxyz(self, name):

        cords = self.cord.frame()
            
        fo = file(name, "w")
        data = self.labels.data
        fo.write("%d\n"%self.natoms())

#   11  HC    -0.023194    1.941312    0.866941     6     9
      
        for i in range(self.natoms()):
            indexer = data[i].field
            fo.write("%6d  %-4.4s%11.6f %11.6f %11.6f %5.5d%s"%(i+1,
                                                                indexer("atomname"),
                                                                cords[i,0],
                                                                cords[i,1],
                                                                cords[i,02],
                                                                indexer("atomtypenum"),
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
            weights = self.labels.data.field("mass")[atomlist]
        # get positions for our atoms
        #   I tested this, it won't alter the original frame.
        frame = self.cord.frame()[atomlist]
        frame.transpose()
        # COM = sum(r*m) / sum(m)
        COM = numarray.sum(frame*weights, axis=1) / sum(weights)

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
