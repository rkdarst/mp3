# Richard Darst, 2005
# pylint: disable-msg=C0101


"""
This module provides functions to get data from PSF files.  

There are two ways to get data out of this file.

- 'bondconnections': This lists what is joined togther NOT DONE
  #FIXME

- 'bonds' or 'angles', etc.

The following options are common to all functions:

- 'number': only return bonds for the first 'number'-1 atoms.  That
    is, the atomnumber you specify is NOT included.  This test occurs
    AFTER 'atomoffset' is processed.  Default to including all atoms.

- 'excludeoutside': If this is is True, then any bonds involving one
    atom excluded are excluded.  The alternative is to just have one
    endpoint returned.

- 'atomnumberoffset': This is is number (integer) to add to each atom
    index.  Use it to compensate for psf indexing from 1, mp3 indexing
    from 0.  (Note that you probably want this to be a negative
    number.)


Dev notes/ To do:
- this function could be made more flexable by returning a
  dictionary instead, and using an exclude list

"""

import itertools

def getnatoms(filename):
    """Return the number of atoms in a PSF

    If '!NATOMS' is not found in the first 100 lines, return None and
    give up.
    """
    i = 0
    for i, line in enumerate(file(filename)):
        if line.find("!NATOM") == -1:
            continue
        if i > 100:
            print "Could not find NATOMS in the first 100 lines, is this a PSF?"
            return None
        natoms = int(line.split()[0])
        return natoms


def getnbonds(filename):
    """Return the number of bonds in a PSF

    If '!NBOND' is not found in the file, return None and
    give up.
    """
    i = 0
    for i, line in enumerate(file(filename)):
        if line.find("!NBOND") == -1:
            continue
        nbonds = int(line.split()[0])
    if nbonds == -1:
        print 'Could not find NBOND in the file.  is this a PSF?'
        return None
    return nbonds


def getnthetas(filename):
    """Return the number of thetas (angles defined by atom triples) in a PSF

    If '!NTHETA' is not found in the file, return None and
    give up.
    """    
    i = 0
    nthetas = -1
    for i, line in enumerate(file(filename)):
        if line.find("!NTHETA") == -1:
            continue
        nthetas = int(line.split()[0])
    if nthetas == -1:
        print 'Could not find NTHETA in the file.  is this a PSF?'
        return None
    return nthetas


def getnphis(filename):
    """Return the number of phis (dihedral angles) in a PSF

    If '!NPHI' is not found in the file, return None and
    give up.
    """    
    i = 0
    nphis = -1
    for i, line in enumerate(file(filename)):
        if line.find("!NPHI") == -1:
            continue
        nphis = int(line.split()[0])
    if nphis == -1:
        print 'Could not find NPHI in the file.  is this a PSF?'
        return None

    return nphis


def getnimphis(filename):
    """Return the number of imphis (improper angles) in a PSF

    If '!NIMPHI' is not found in the file, return None and
    give up.
    """    
    i = 0
    nimphis = -1
    for i, line in enumerate(file(filename)):
        if line.find("!NIMPHI") == -1:
            continue
        nimphis = int(line.split()[0])
    if nimphis == -1:
        print 'Could not find NIMPHI in the file.  is this a PSF?'
        return None

    return nimphis


def getatomlist(filename, nprotatoms=0):
    """Returns a list of dictionaries of all the atoms in the .psf file and containing all
    of the .psf information in the atom section.  Can write to the Atom Defs section of
    a .top file"""
    atomlist = []
    inFile = file(filename)
    for line in inFile:
        if line.find('!NATOM') == -1:
            continue
        if nprotatoms == 0:
            nprotatoms = int(line.split()[0])        
        break

    i = 0
    for line in inFile:
        if line.strip() == "":
            break
        if i >= nprotatoms:
            break
        atomlist.append({})
        atomindex = line.split()[0]
        resgroup = line.split()[3]
        atomtype = line.split()[5]
        charge = line.split()[6]
        mass = line.split()[7]
        atomlist[i]['atomindex'] = atomindex
        atomlist[i]['resgroup'] = resgroup
        atomlist[i]['atomtype'] = atomtype
        atomlist[i]['charge'] = charge
        atomlist[i]['mass'] = mass
        i += 1


    return atomlist


def getbondlist(filename):
    """Returns a tuple that contains the number of bonds as one
    element of the tuple and a list of bonded pairs of atoms as the
    other element."""
    
    bondlist = []        
    inFile = file(filename)
    for line in inFile:
        if line.find('!NBOND') == -1:
            continue
        nbonds = int(line.split()[0])
        break

    # the .psf file lists four pairs of atoms per line
    for line in inFile:
        if line.strip() == "":
            break
        for i in xrange(len(line.split())/2):
            bondlist.append((line.split()[2*i],line.split()[2*i+1]))
    
    return nbonds,bondlist



def getbondconnections(filename, number=None, excludeoutside=True, atomnumberoffset=0):
    """Return the list of connections via bonds in a protein.
    """
    # get number of atoms, set number of bonds to return
    natoms = getnatoms(filename)
    if number is None:
        number = natoms
    bonds = [ [] for i in xrange(number) ]
    fo = file(filename)
    # skip natoms
    for i in xrange(natoms):
        fo.readline(50)
    # find nbonds
    for line in fo:
        if line.find("!NBOND") == -1:
            continue
        nbonds = int(line.split()[0])
        break
    # parse the bonds:
    allbonds = [ ]
    for line in fo:
        #print line
        if line.strip() == "":
            break
        allbonds.extend( [int(i)+atomnumberoffset for i in line.split() ] )
        if excludeoutside:
        # skip it if either is outside the range
            for i in range(len(allbonds)/2):
                a, b = allbonds[0:2]
                del allbonds[0:2]
                if a < number and b < number:
                    #continue
                    bonds[a].append(b)
                    bonds[b].append(a)
        else:
        # include just one endpoint in the range
            for i in range(len(allbonds)/2):
                a, b = allbonds[0:2]
                del allbonds[0:2]
                if a < number:
                    bonds[a].append(b)
                if b < number:
                    bonds[b].append(a)

    return {"n": nbonds,
            "connections": bonds }


def getthetalist(filename):
    """Returns a list of triples defining and theta.
    """
    thetalist = []        
    inFile = file(filename)
    for line in inFile:
        if line.find('!NTHETA') == -1:
            continue
        ntheta = int(line.split()[0])
        break

    for line in inFile:
        if line.strip() == "":
            break
        
        # the .psf file lists three triples of atoms per line        
        for i in xrange(len(line.split())/3):      
            thetalist.append((line.split()[3*i],line.split()[3*i+1],line.split()[3*i+2]))

    inFile.close()
    return ntheta,thetalist



def getangleconnections(filename, number=None, excludeoutside=True, atomnumberoffset=0):
    """Return the connections via angles in a PSF.
    """
    # get number of atoms, set number of bonds to return
    natoms = getnatoms(filename)
    if number is None:
        number = natoms
    connections = [ [] for i in xrange(number) ]
    fo = file(filename)
    # skip natoms
    for i in xrange(natoms):
        fo.readline(50)
    # find nbonds
    for line in fo:
        if line.find("!NTHETA") == -1:
            continue
        nangles = int(line.split()[0])
        break
    # parse the bonds:
    allbonds = [ ]
    for line in fo:
        #print line
        if line.strip() == "":
            break
        allbonds.extend( [int(i)+atomnumberoffset for i in line.split() ] )
        if excludeoutside:
        # skip it if either is outside the range
            for i in range(len(allbonds)/3):
                a, b, c = allbonds[0:3]
                del allbonds[0:3]
                if a < number and b < number and c < number:
                    #continue
                    connections[a].append(b)
                    connections[a].append(c)
                    connections[b].append(a)
                    connections[b].append(c)
                    connections[c].append(a)
                    connections[c].append(b)
        else:
        # include just one endpoint in the range
            for i in range(len(allbonds)/2):
                a, b = allbonds[0:2]
                del allbonds[0:2]
                if a < number:
                    connections[a].append(b)
                    connections[a].append(c)
                if b < number:
                    connections[b].append(a)
                    connections[b].append(c)
                if c < number:
                    connections[c].append(a)
                    connections[c].append(b)

    return {"n": nangles,
            "connections": connections }


def getdihedrallist(filename):
    """Returns a list of dihedral forming atoms."""
    dihedrallist = []        
    inFile = file(filename)
    for line in inFile:
        if line.find('!NPHI') == -1:
            continue
        nphi = int(line.split()[0])
        break

    for line in inFile:
        if line.strip() == "":
            break
        
        # the .psf file lists two quadruples of atoms per line        
        for i in xrange(len(line.split())/4):      
            dihedrallist.append((line.split()[4*i],line.split()[4*i+1],line.split()[4*i+2],
                                 line.split()[4*i+3]))

    inFile.close()
    return nphi,dihedrallist



def getdihedralconnections(filename, number=None, excludeoutside=True, atomnumberoffset=0):
    """Return the connections via dihedrals in a PSF.
    """
    # get number of atoms, set number of bonds to return
    natoms = getnatoms(filename)
    if number is None:
        number = natoms
    connections = [ [] for i in xrange(number) ]
    fo = file(filename)
    # skip natoms
    for i in xrange(natoms):
        fo.readline(50)
    # find nbonds
    for line in fo:
        if line.find("!NPHI") == -1:
            continue
        ndihedrals = int(line.split()[0])
        break
    # parse the bonds:
    allbonds = [ ]
    for line in fo:
        #print line
        if line.strip() == "":
            break
        allbonds.extend( [int(i)+atomnumberoffset for i in line.split() ] )
        if excludeoutside:
        # skip it if either is outside the range
            for i in range(len(allbonds)/4):
                a, b, c, d = allbonds[0:4]
                del allbonds[0:4]
                if a < number and b < number and c < number and d < number:
                    #continue
                    connections[a].append(b)
                    connections[a].append(c)
                    connections[a].append(d)

                    connections[b].append(a)
                    connections[b].append(c)
                    connections[b].append(d)
                    
                    connections[c].append(a)
                    connections[c].append(b)
                    connections[c].append(d)
                    
                    connections[d].append(a)
                    connections[d].append(b)
                    connections[d].append(c)
        else:
        # include just one endpoint in the range
            for i in range(len(allbonds)/2):
                a, b = allbonds[0:2]
                del allbonds[0:2]
                if a < number:
                    connections[a].append(b)
                    connections[a].append(c)
                    connections[a].append(d)
                if b < number:
                    connections[b].append(a)
                    connections[b].append(c)
                    connections[b].append(d)
                if c < number:
                    connections[c].append(a)
                    connections[c].append(b)
                    connections[c].append(d)
                if d < number:
                    connections[d].append(a)
                    connections[d].append(b)
                    connections[d].append(c)

    return {"n": ndihedrals,
            "connections": connections }


def getimproperlist(filename):
    """Returns a list of improper forming atoms."""
    improperlist = []        
    inFile = file(filename)
    for line in inFile:
        if line.find('!NIMPHI') == -1:
            continue
        nimphi = int(line.split()[0])
        break

    for line in inFile:
        if line.strip() == "":
            break
        
        # the .psf file lists two quadruples of atoms per line        
        for i in xrange(len(line.split())/4):      
            improperlist.append((line.split()[4*i],line.split()[4*i+1],line.split()[4*i+2],
                                 line.split()[4*i+3]))

    inFile.close()
    return nimphi,improperlist




def getimproperconnections(filename, number=None, excludeoutside=True, atomnumberoffset=0):
    """Return the connections via impropers in a PSF.
    """
    # get number of atoms, set number of bonds to return
    natoms = getnatoms(filename)
    if number is None:
        number = natoms
    connections = [ [] for i in xrange(number) ]
    fo = file(filename)
    # skip natoms
    for i in xrange(natoms):
        fo.readline(50)
    # find nbonds
    for line in fo:
        if line.find("!NIMPHI") == -1:
            continue
        ndihedrals = int(line.split()[0])
        break
    # parse the bonds:
    allbonds = [ ]
    for line in fo:
        #print line
        if line.strip() == "":
            break
        allbonds.extend( [int(i)+atomnumberoffset for i in line.split() ] )
        if excludeoutside:
        # skip it if either is outside the range
            for i in range(len(allbonds)/4):
                a, b, c, d = allbonds[0:4]
                del allbonds[0:4]
                if a < number and b < number and c < number and d < number:
                    #continue
                    connections[a].append(b)
                    connections[a].append(c)
                    connections[a].append(d)

                    connections[b].append(a)
                    connections[b].append(c)
                    connections[b].append(d)
                    
                    connections[c].append(a)
                    connections[c].append(b)
                    connections[c].append(d)
                    
                    connections[d].append(a)
                    connections[d].append(b)
                    connections[d].append(c)
        else:
        # include just one endpoint in the range
            for i in range(len(allbonds)/2):
                a, b = allbonds[0:2]
                del allbonds[0:2]
                if a < number:
                    connections[a].append(b)
                    connections[a].append(c)
                    connections[a].append(d)
                if b < number:
                    connections[b].append(a)
                    connections[b].append(c)
                    connections[b].append(d)
                if c < number:
                    connections[c].append(a)
                    connections[c].append(b)
                    connections[c].append(d)
                if d < number:
                    connections[d].append(a)
                    connections[d].append(b)
                    connections[d].append(c)

    return {"n": ndihedrals,
            "connections": connections }



def commonconnections(*lists):
    """'zips' a list of connections.

    'zips' a list of connections returned by other connection
    functions above.

    NOTE: they must have the same indexes...
    """

    maxlen = max( [len(list_) for list_ in lists] )
    connections = [ [] for i in xrange(maxlen) ]
    for list_, otherlists in itertools.izip(connections, itertools.izip(*lists)):
        for atom in itertools.chain(*otherlists):
            atom in list_ or list_.append(atom)

    return connections



if __name__ == "__main__":
    if False:
        filename = "files/apomb2C_SPCEions.psf"
        print getnatoms(filename)
        bondconn = getbondconnections(filename,
                                 number=2464,
                                 excludeoutside=True,
                                 atomnumberoffset=-1)
        angleconn = getangleconnections(filename,
                                  number=2464,
                                  excludeoutside=True,
                                  atomnumberoffset=-1)
        dihedralconn = getdihedralconnections(filename,
                                     number=2464,
                                     excludeoutside=True,
                                     atomnumberoffset=-1)
        #print bondconn
        #print angleconn
        #print dihedralconn

        print commonconnections(bondconn["connections"],
                                angleconn["connections"],
                                dihedralconn["connections"])

    if True:
##         Tests added by John R. Dowdle, 04/24/2007 to test the
##         getatomlist, getbondlist, getthetalist, getdihedrallist,
##         and getimproperlist functions.
        filename = 'files/ubq_patched.psf'
        atomlist = getatomlist(filename,nprotatoms=1231)
        print atomlist
        nbonds,bondlist = getbondlist(filename)
        print 'nbonds = ',nbonds
        print 'len(bondlist) = ',len(bondlist)
        print 'first five bonds:  '
        print bondlist[0:5]
        print 'last five bonds:  '
        print bondlist[len(bondlist)-5:len(bondlist)]        
        
        nthetas,thetalist = getthetalist(filename)
        print 'nthetas = ',nthetas
        print 'len(thetalist) = ',len(thetalist)
        print 'first five thetas:  '
        print thetalist[0:5]
        print 'last five thetas:  '
        print thetalist[len(thetalist)-5:len(thetalist)]        
        
        ndihedrals,dihedrallist = getdihedrallist(filename)
        print 'ndihedrals = ',ndihedrals
        print 'len(dihedrallist) = ',len(dihedrallist)
        print 'first five dihedrals:  '
        print dihedrallist[0:5]
        print 'last five dihedrals:  '
        print dihedrallist[len(dihedrallist)-5:len(dihedrallist)]        
        
        nimpropers,improperlist = getimproperlist(filename)
        print 'nimpropers = ',nimpropers
        print 'len(improperlist) = ',len(improperlist)
        print 'first five impropers:  '
        print improperlist[0:5]
        print 'last five impropers:  '
        print improperlist[len(improperlist)-5:len(improperlist)]        
        
        nbonds = getnbonds(filename)
        print 'nbonds = ',nbonds
        nthetas = getnthetas(filename)
        print 'nthetas = ',nthetas
        nphis = getnphis(filename)
        print 'nphis = ',nphis
        nimphis = getnimphis(filename)
        print 'nimphis = ',nimphis
        
