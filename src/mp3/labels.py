import string
import glob
import sets
import numarray.records
import logging
import os.path
import mp3.log
import mp3.cm3dparse as cm3dparse
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

    examples:
        list with masses of all atoms:
        L.data.field("mass")


    """

    def natoms(self):
        return self._natoms

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

        Loads metadata from only the first PDB given.

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
            natoms = self._natoms

        lists_to_intersect = []
        for key,value in keywords.iteritems():
            # We want to be able to pick if we are finding an exact
            # match, or a range.  If it is an exact match, then
            # "value" will be the element to match.  If it is a range,
            # "value" should be a tuple containing (low, high),
            # inclusively.

            fiel = self.data.field(key)
            if type(value) == tuple and len(value) == 2:  # we have a range match
                low = value[0]
                high = value[1]
                #new_list = [ i for i in range(natoms)
                #             if ( self.data.field(key)[i] >= low and self.data.field(key)[i] <= high ) ]
                new_list = [ i for i in range(natoms)
                             if ( fiel[i] >= low and fiel[i] <= high ) ]
                lists_to_intersect.append(new_list)
            else:                     # we have an exact match
                new_list = [ i for i in range(natoms)
                              if fiel[i] == value ]
            # End branch between exact and exact matches.
            lists_to_intersect.append(new_list)
        # this should be moved the built in `set`, implemented in C.
        # But it's only in python2.4.
        atomset = sets.Set(lists_to_intersect[0] )  
        for otherlist in lists_to_intersect[1:]:
            atomset &= sets.Set(otherlist)
        return atomset
        #return atomlist

    findrange = findatoms

    def _makedataarray(self):
        """Makes the numarray record array to store labels in.
        """
        
        self.data = numarray.records.array(formats='a4,Int32,a4,a4,a4,Float32,Float32,Int32,Int32,Int32,a2,Int32', names='segname,resnum,resname,atomname,atomtype,charge,mass,unused,occupancy,tempfactor,element,atomtypenum', shape=self._natoms, buffer='\000'*(11*4+2)*self._natoms)


    def _parsepsf(self, psffileobject):
        """Loads a psf from the given file-object.

        Not a public method, check out getfrompsf() instead.
        """
        if hasattr(self, "data") == True:
            logging.error("We already have data ...overwriting with the new stuff.")
    
    
        line = psffileobject.readline()              #read and verify the "psf" on the first line
        if line[0:3] != "PSF":
            logging.error("Is this a psf file?\n")
    
            #print "hi"
                
        line = psffileobject.readline()              #find the NATOMS line and read it in
        while line.find("NATOM") == -1:
            line = psffileobject.readline()
        self._natoms = int(line.split()[0])
    
        if not hasattr(self, "data"):
            self._makedataarray()
        
        for atomn in range(0,self._natoms):
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
        natoms = 0
        while True:
            line = pdbfo.readline()
            if line[0:4] == "ATOM":
                break
        natoms += 1
        while True:
            line = pdbfo.readline()
            if line[0:4] == "ATOM":
                natoms += 1
            else:
                break
            
        self._natoms = natoms
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
            if atomn+1 > self._natoms:  # bug city!
                         #let's be sure we get this right. when atomn+1 = natoms,
                         # we are where we need to be to stop. But we still have to
                         # go around one more time to be sure to get the current atom.
                break
            line = pdbfo.readline()


    def _parsetinkerxyz(self, xyzfo):
        """Reads atomname and atomtype from a tinker .xyz file object
        """
        self._natoms = int(xyzfo.readline().strip() )
        if not hasattr(self, "data"):
            self._makedataarray()
        self._makebondlist()
        for atomn in xrange(self._natoms):
            line = xyzfo.readline()
            self.data.field('atomname')[atomn] = line[7:11].strip()
            self.data.field('atomtypenum')[atomn] = int(line[48:53])
            self._tinkerbondlist.append(line[53:])

    def _makebondlist(self):
        self._tinkerbondlist =  []


    def readcm3dset(self, setfile,
                    fields=('atomtype', 'mass', 'charge', 'resname','segname')):
        """Read in labels information from a CM3D set file.

        `setfile` -- set file to read.

        `fields` -- optional argument telling which fields we want to
        store/override in the labels data.  This should be a tuple
        with any combination of the strings 'atomtype', 'mass',
        'charge', 'resname', 'segname', or 'atomname' in it.
        Specifing less will slightly speed it up.

        Mapping of CM3D fields --> labels fields
        mol_index in setfile   --> segname
        atom_typ  in topfile   --> atomtype
        mass      in topfile   --> mass
        charge    in topfile   --> change
        group     in topfile   --> group
        sttype    in topfile       (not used)
        atomname  in topfile * --> atomname

        * `atomname` is a non-standard field.  If you code it with
          some special syntax, mess with the regular expression in
          cm3dparse.rKW to make it detect it.
        """
        mp3.log.info("in Labels.readcm3dset, loading %s"%setfile)
        dirname = os.path.dirname(setfile)
        natoms = 0

        # Preselect and store the compiling command
        codeLine = ""
        if 'atomtype' in fields:
            codeLine += 'data.field("atomtype")[N] = contents_["atom_typ"].strip() \n'
        if 'mass' in fields:
            codeLine += 'data.field("mass")[N] = float(contents_.get("mass", 0.)) \n'
        if 'charge' in fields:
            codeLine += 'data.field("charge")[N] = float(contents_.get("charge", 0.)) \n'
        if 'resname' in fields:
            codeLine += 'data.field("resname")[N] = contents_.get("group", "").strip() \n'
        if 'segname' in fields:
            # MOL_INDEX is a integer, but we store it as a string in the segname field.
            codeLine += 'data.field("segname")[N] = MOL_INDEX \n'
        if 'atomname' in fields:
            codeLine += 'data.field("atomname")[N] = contents_.get("atomname", "").strip() \n'
        # contents_["sttype"] is unused
        print 'the code line'
        codeLine = compile(codeLine,'<custom_command>','exec')
        
        # First, read through all files once to figure out how many atoms there are.
        setfiledata = file(setfile).read()
        for MKW, contents in cm3dparse.iterMKWd(setfiledata):
            if MKW != "mol_def":
                continue
            nmol = int(contents["nmol"])  # number of times that molecule is repeated
            molParamFile = os.path.join(dirname, contents["mol_parm_file"])
            topfiledata = file(molParamFile).read()
            molNAtom = None # for debugging purposes below
            for MKW_, contents_ in cm3dparse.iterMKWd(topfiledata):
                if MKW_ != "mol_name_def":
                    continue
                molNAtom = int(contents_["natom"])
                break
            if molNAtom is None:
                mp3.log.error("We didn't find mol_name_def in the topology file %s"%molParamFile)
            mp3.log.info("Loading file %s, have %s replicates of %s atoms"%(molParamFile,
                                                                            nmol,
                                                                            molNAtom))
            natoms += nmol * molNAtom
        mp3.log.info("We have found %s atoms total in %s"%(natoms, setfile))
        self._natoms = natoms

        # Begin loading in all the other data.
        if not hasattr(self, "data"):
            self._makedataarray()
        data = self.data
        # Now, we have to iterate through all of that again.  I don't
        # like to just copy the code from above, but I guess I'll do
        # that.  We can streamline it since presumably it's
        # well-behaved (has the required fields) now.
        N = 0   # atom index we are currently on.
        for MKW, contents in cm3dparse.iterMKWd(setfiledata):
            #print MKW, contents
            if MKW != "mol_def":
                continue
            MOL_INDEX = contents["mol_index"]
            molParamFile = os.path.join(dirname, contents["mol_parm_file"])
            topfiledata = file(molParamFile).read()
            nmol = int(contents["nmol"])  # number of times that molecule is repeated
            for i in range(nmol):
                for MKW_, contents_ in cm3dparse.iterMKWd(topfiledata):
                    #print MKW_, contents_
                    if MKW_ == "mol_name_def":
                        molNAtom = int(contents_["natom"])
                        #ResName = contents_["mol_name"]
                        MOL_NAME = contents_["mol_name"]
                    if MKW_ == "atom_def":
                        exec codeLine
                        #print N
                        molNAtom -= 1
                        N += 1
                    if molNAtom <= 0:
                        break
