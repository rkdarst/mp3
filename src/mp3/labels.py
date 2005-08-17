import string
import glob
import sets
import numarray.records
import logging
thelog = logging.getLogger('mp3')

class Labels:
    """This class handles extra atom information, such as names and residue numbers.

    It reads information from psf files. Note that it isn't very smart at
    parsing these psf's yet.

    Class contians:
        getfrompsf() -- Gets label information from a psf
        getfrompdb() -- Gets label information from a pdb
        getfromtxyz() -- Gets labess from tinker xyz
        findatoms() -- Search for atoms matching given criteria
        findrange() -- Same as findrange
        data -- Numarray record array containing actual labels.  It is
                indexed by atomnumber, starting from zero
            field() -- method which accesses fields of the data array.
                field names (atomnum stored in field index, all other
                             fields correspond to the values found in
                             the PDB/PSF/etc. ):
                    atomtype -- string4
                    atomtypenum -- int   (from tinker)
                    atomname -- string4
                    charge -- float
                    element -- string2
                    mass -- float
                    occupancy -- float
                    resname -- string4
                    resnum  -- int
                    segname -- string4
                    tempfactor -- float
                    unused  -- int
      
    """

    def getfrompsf(self, psfname):
        """Loads atom information from the given psf file.

        Pass this function a filename, and it will load labels from it.
        """
        #should open the file if it is a string..., otherwise pass the file object through
        psffileobject = file(psfname,"r")
        self._parsepsf(psffileobject)
        psffileobject.close()

    def getfrompdb(self, pdblist):
        """Loads labels from a PDB.

        You can pass this function a string or a list-like object. If it
        is a string, it is expanded by shell-glob rules and sorted, and
        the first element is used.  Otherwise, it must be a list-like
        object and the first element must be a string and it is used.
        """

        if type(pdblist) == str:
            pdblist = glob.glob(pdblist)
            pdblist.sort()
        self._parsepdb(file(pdblist[0],"r"))

    def getfromtxyz(self, tinkerxyzlist):
        """Load labels from a Tinker XYZ file.

        This operates analogously to the getfrompdb method.
        """
        if type(tinkerxyzlist) == str:
            tinkerxyzlist = glob.glob(tinkerxyzlist)
            tinkerxyzlist.sort()
        self._parsetinkerxyz(file(tinkerxyzlist[0],"r"))

    def findatoms(self, **keywords):
        """This function finds atoms matching specified conditions.

        This method takes keyword argument(s) of the form FIELD=VALUE.
        FIELD is a fieldname of our data array, as given in the class
        docstring (help(mp3.Labels)).  

        VALUE can be an exact match, or a range of matches.  If VALUE
        is a tuple of length 2, it is interpreted as the range
        VALUE=(low, high).  The criteria is matched if
        "low <= field <= high".  If VALUE is a single object (such as
        5 or"H2"), the criteria is match if "field == VALUE".

        A "set" of atoms is returned which contains all atoms matching all
        criteria (an AND operator).

        The atoms are returned in a "set" object.  See documentation on
        the standard python module "sets" to find properties of these
        objects.  It ensures that there are no duplicates, but does not
        preserve order.  To convert this into a list, use list(your_set).
        If the list must be sorted, used "thelist.sort()" (which
        operates in-place, so do "thelist.sort()", not
        "thelist=thelist.sort()").

        If there is an "natoms=INTEGER" keyword, restrict search to
        the first INTEGER atoms.  This can be useful when, for
        example, you are searching for a protein which is the first
        3000 atoms of a 100000 atom file, and you want to only search
        the first 3000 atoms, not all 100000 (for time reasons).
        """
        #atomlist = []
        #for i in range(self.natoms):
        #    if self.data.field(field)[i] == name:
        #        atomlist.append(i)

        # This code is to advance to a new, set-based method of
        # selecting atoms.  I haven't switched to it yet.
        # This should take an intersection of all given conditions.  In
        # order to do this, I will find a list of atoms which match each
        # condition, and then find the intersection of all of them.

        # Process natoms: if not given, natoms=self.natoms(), if given,
        # only iterate over those many atoms

        if keywords.has_key("natoms"):
            natoms = keywords["natoms"]
            del keywords["natoms"]
        else:
            natoms = self.natoms

        lists_to_intersect = []
        for key,value in keywords.iteritems():
            # We want to be able to pick if we are finding an exact
            # match, or a range.  If it is an exact match, then
            # "value" will be the element to match.  If it is a range,
            # "value" should be a tuple containing (low, high),
            # inclusively.

            if type(value) == tuple and len(value) == 2:  # we have a range match
                low = value[0]
                high = value[1]
                fiel = self.data.field(key)
                #new_list = [ i for i in range(natoms)
                #             if ( self.data.field(key)[i] >= low and self.data.field(key)[i] <= high ) ]
                new_list = [ i for i in range(natoms)
                             if ( fiel[i] >= low and fiel[i] <= high ) ]
                lists_to_intersect.append(new_list)
            else:                     # we have an exact match
                new_list = [ i for i in range(natoms)
                              if self.data.field(key)[i] == value ]
            # End branch between exact and exact matches.
            lists_to_intersect.append(new_list)
        atomset = sets.Set(lists_to_intersect[0] )
        for otherlist in lists_to_intersect[1:]:
            atomset &= sets.Set(otherlist)
        return atomset
        #return atomlist

    findrange = findatoms

    def _makedataarray(self):
        """Makes the numarray record array to store labels in.
        """
        
        self.data = numarray.records.array(formats='a4,Int32,a4,a4,a4,Float32,Float32,Int32,Int32,Int32,a2,Int32', names='segname,resnum,resname,atomname,atomtype,charge,mass,unused,occupancy,tempfactor,element,atomtypenum', shape=self.natoms, buffer='\000'*(11*4+2)*self.natoms)


    def _parsepsf(self, psffileobject):
        """Loads a psf from the given file-object.

        Not a public method, check out getfrompsf() instead.
        """
        if hasattr(self, "data") == True:
            logging.error("We already have data! ...overwriting")
    
    
        line = psffileobject.readline()              #read and verify the "psf" on the first line
        if line[0:3] != "PSF":
            logging.error("Is this a psf file?\n")
    
            #print "hi"
                
        line = psffileobject.readline()              #find the NATOMS line and read it in
        while line.find("NATOM") == -1:
            line = psffileobject.readline()
        self.natoms = int(line.split()[0])
    

        self._makedataarray()
        
        for atomn in range(0,self.natoms):
            #print atomn
            line = psffileobject.readline()
            #fields= line.split()
    
            #if int(fields[0]) -1 != atomn:
                #print natom 
            #    thelog.error("Atom number mismatch at %d\n" % (atomn + 1) )
              
    
            #self.data.field('segname')[atomn] = fields[1].strip()
            #self.data.field('resnum')[atomn] = int(fields[2])
            #self.data.field('resname')[atomn] = fields[3].strip()
            #self.data.field('atomname')[atomn] = fields[4].strip()    #used to be "atomid"
            #self.data.field('atomtype')[atomn] = fields[5].strip()
            #self.data.field('charge')[atomn] = float(fields[6])
            #self.data.field('mass')[atomn] = float(fields[7])
            #self.data.field('unused')[atomn] = int(fields[8])  #used to be "junk"

#            print self.data.field('atomname')[0]

            self.data.field('segname')[atomn] = line[9:13].strip()
            self.data.field('resnum')[atomn] = int(line[14:18])
            self.data.field('resname')[atomn] = line[19:23].strip()
            self.data.field('atomname')[atomn] = line[24:28] #not stripped!! ... spacing must be preserved    #used to be "atomid"
            self.data.field('atomtype')[atomn] = line[29:33].strip()
            self.data.field('charge')[atomn] = float(line[34:44])
            self.data.field('mass')[atomn] = float(line[44:58])  #not sure about this spacing
            self.data.field('unused')[atomn] = int(line[69])  #used to be "junk"
    
            #print fields[1]



    def _parsepdb(self, pdbfo):
        """Reads labels from a PDB file-object

        Use _parsepdb instead.  The argument is a open file-object.
        """

        ### we've got to find out how many atoms there are.
          #  it would be better if we could get the data from the
          #  cord object, but no one says that they will be related.
          #  it's more flexable (and you only have to do it once) to
          #  do it here independently
          # code copied straight from pdbcord.py
        self.natoms = 0
        while True:
            line = pdbfo.readline()
            if line[0:4] == "ATOM":
                break
        self.natoms += 1
        while True:
            line = pdbfo.readline()
            if line[0:4] == "ATOM":
                self.natoms += 1
            else:
                break

        self._makedataarray()
        ### now that we know the number of atoms, go parse the thing
        pdbfo.seek(0)
        while True:
            line = pdbfo.readline()
            if line[0:4] == "ATOM":   # we find the first atom line here, but
                                      #don't process it until the next loop.
                break
        atomn = 0
        while True:
            self.data.field('atomname')[atomn] = line[12:16]
            self.data.field('resname')[atomn] = line[17:21].strip()  # illegally expand this field
                                                                     # by one
            self.data.field('resnum')[atomn] = int(line[22:26])
            self.data.field('occupancy')[atomn] = float(line[54:60])
            self.data.field('tempfactor')[atomn] = float(line[60:66])
            self.data.field('segname')[atomn] = line[72:76].strip()
            self.data.field('element')[atomn] = line[76:78]  #used to be "junk"


            atomn += 1
            if atomn+1 > self.natoms:  # bug city!
                         #let's be sure we get this right. when atomn+1 = natoms,
                         # we are where we need to be to stop. But we still have to
                         # go around one more time to be sure to get the current atom.
                break
            line = pdbfo.readline()


    def _parsetinkerxyz(self, xyzfo):
        """Reads atomname and atomtype from a tinker .xyz file object
        """
        self.natoms = int(xyzfo.readline().strip() )
        self._makedataarray()
        self._makebondlist()
        for atomn in xrange(self.natoms):
            line = xyzfo.readline()
            self.data.field('atomname')[atomn] = line[7:11].strip()
            self.data.field('atomtypenum')[atomn] = int(line[48:53])
            self._tinkerbondlist.append(line[53:])

    def _makebondlist(self):
        self._tinkerbondlist =  []

