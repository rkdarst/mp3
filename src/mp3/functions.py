# This module contains functions to simplify various processing
# using mp3.
#
#
#

import numarray, math
import re
import mp3
import mp3.log

def pdbsystem(name):
    """Easy way to read in PDB's cord and labels.

    The argument is a string or list (same rules as the other pdb
    functions) and returns a <system> object with labels read and
    cord object initilized and moved to the first (perhaps only)
    frame.
    """
    system = mp3.System()
    system.labels.getfrompdb(name)
    system.cord = mp3.CordPDB()
    system.cord.pdblist(name)
    system.cord.nextframe()
    return system


#
#  The "smartsuite"
#
#  all of these are designed to ease the opening of files.  


def smartsystem(*args):
    """Create a system from arbitrary input files.

    This function takes any arguments, and returns a system.  The
    systems contains all coordinate files, in sequence, as well as any
    label files.
    """
    inputs = []
    # Make a single list out of all inputs, regardless if the
    # arguments are lists or single strings.
    for i in args:  
        if type(i) == list or type(i) == tuple:
            inputs.extend(i)
        else:
            inputs.append(i)

    cords = [ ]
    labels = [ ]
    alreadyfoundlabels = []   # we don't want to get labels from multiple tinker files
    for input_ in inputs:
        type_ = whatisit(input_)
        if type_ in ("tinkerxyz", "pdb", "dcd", "tinkerarc"):
            cords.append(input_)
        if type_ == "psf":
            labels.append(input_)
        if type_ in ("tinkerxyz", "tinkerarc", "pbd") and not type_ in alreadyfoundlabels:
            alreadyfoundlabels.append(type_)
            labels.append(input_)
    mp3.log.info("smartsystem: found %d cords, %d labels."%
                 (len(cords), len(labels)))
        
    S = mp3.System()
    if len(cords) == 0:                      # Set up all the cords
        pass
    else:
        C = smartcord(cords)
        S.setcord(C)
    for label in labels:                   # Set up all the labels
        type_ = whatisit(label)
        if type_ == "psf":
            S.labels.getfrompsf(label)
        if type_ == "pdb":
            S.labels.getfrompdb(label)
        if type_ in ("tinkerxyz", "tinkerarc"):
            S.labels.getfromxyz(label)
    #if len(labels) > 0:
    #    print "we don't support reading in labels yet! (pester the maintainer)"

    return S
        

def smartcord(*args):
    """Return Turn the list of inputs into a cord object.
    """
    inputs = []
    for i in args:                               # make it all be one big list of strings
        if type(i) == list or type(i) == tuple:
            inputs.extend(i)
        else:
            inputs.append(i)
    cords = []            # cords is a list of:
                           # strings
                           # lists of strings of pdbs or tinkerxyzs.
    lasttype = None
    for input_ in inputs: # now, make the xys's and pdb's be lists of
                                                 # themselves
        type_ = whatisit(input_)
        if type_ in ("tinkerxyz", "pdb"):
            if lasttype == type_:
                cords[-1].append(input_)
            else:
                cords.append( [ input_ ] ) # append a LIST here, since
                                            # it might be appended to
                                            # later
            lasttype = type_
        elif type_ in ("dcd", "tinkerarc"):
            cords.append(input_)
            lasttype = None
        else:
            print "Unknown file type (skipping):", input_
    #print cords
    # we have better not have a null input list!
    if len(cords) == 0:
        # this is bad...
        return None
    elif len(cords) == 1:
        C = _smartcord_single(cords[0])
    else:
        cord_objs = [ _smartcord_single(i) for i in cords ]
        print cord_objs
        C = mp3.CordMerge(cords=cord_objs)
    return C

def _smartcord_single(name):
    """Return a cord object.

    Takes either a string with a cord file, or a list of tinkerxyz's
    or pdb's.
    """
    type_ = whatisit(name)
    #print name, "lalala"
    if type_ in ("tinkerxyz", "tinkerxyzlist"):

        # if the length of the list is one, what if it is really a
        # glob (like "*.xyz_*") that the user was expecting the
        # CordTinkerXYZ to expand as a shell glob.  It will only do
        # this if you pass it a string, not a list with stuff in it.
        #
        # Note that there is still one case where this can break.  If
        # the user passed two globs in a row, then it won't turn them
        # into two individual strings, it would leave it as a list
        # with two (non-globbed) strings in it.  Let the user beware
        # (not that I'm warning them any, yet).
        if len(name) == 1:
            return mp3.CordTinkerXYZ(xyzlist=name[0])
        else:
            return mp3.CordTinkerXYZ(xyzlist=name)
    if type_ in ("tinkerarc",):
        return mp3.CordTinkerArc(arc=name)
    if type_ in ("pdb", "pdblist"):
        if len(name) == 1:
            return mp3.CordPDB(pdblist=name[0])
        else:
            return mp3.CordPDB(pdblist=name)
    if type_ in ("dcd",):
        return mp3.CordDCD(dcd=name)

# filename_regex is a mapping between file types and (compiled)
# regular expressions on the file name that tell what kind of files
# they are.  You can add keys to it if you want.
# you should SEARCH (not match) for these rexexps

filename_regex = {}
filename_regex["tinkerxyz"] = re.compile(r'\.xyz(_\d+)?$')
filename_regex["tinkerarc"] = re.compile(r'\.arc(_\d+)?$')
filename_regex["dcd"]       = re.compile(r'\.dcd$')
filename_regex["pdb"]       = re.compile(r'\.pdb$')
filename_regex["psf"]       = re.compile(r'\.psf$')

#mp3.functions.filename_regex["tinkerxyz"] = re.compile(r'\.(\d\d\d+)$')


def whatisit(filename):
    """Return the type of file.

    Do a simple regex match on the given filename.  Return a string
    indicating what type of file it probably is.  Return values are
    "tinkerxyz", "tinkerarc", "dcd", "pdb", "psf", or None (if it
    doesn't match any of them).

    Can also return "pdblist" or "tinkerxyzlist".  Does NOT handle
    non-homogeneous or zero-length lists.
    """
    # you should SEARCH (not match) for these rexexps
    regex = filename_regex
    for key, value in regex.iteritems():
        if type(filename) == list:
            if regex[key].search(filename[0]) != None:
                return key+"list"
        else:
            if regex[key].search(filename) != None:
                #print key
                return key
            
    return None
    
    
def rmsd(frame1, frame2):
    """Finds the RMSD between two frames

    """
    ## find the displacement for each coordinate
    disp = frame1 - frame2
    ##find squared displacement
    numarray.multiply(disp, disp, disp)
    sd = numarray.sum(disp, axis=1)    # "squared displacement"
    ## mean squared displacement, etc.
    msd = sd.sum() / len(sd)
    rmsd = math.sqrt(msd)

#    print msd
#    print rmsd
    return rmsd


def cordtransform(frame, move=None, rotate=None):
    """transforms a frame in place.

    This function takes two three-vectors, given by keyword. Angles,
    of course, are in radians.

    ! The frame is rotated before it is translated.  Rotations are
    done in the order x,y,z !

    move=(x-translate, y-translate, z-translate)
    rotate=(x-rotate, y-rotate, z-rotate)

    Returns the new frame.  

    """
    sin = math.sin
    cos = math.cos

    ## Rotate the frame.
    if rotate != None:
        ## first task: create our net rotation matrix:
        xmatrix = [ [ 1,                0,              0  ],
                    [ 0,   cos(rotate[0]),  sin(rotate[0]) ],
                    [ 0,  -sin(rotate[0]),  cos(rotate[0]) ], ]

        ymatrix = [ [ cos(rotate[1]),  0,  -sin(rotate[1]) ],
                    [              0,  1,               0  ],
                    [ sin(rotate[1]),  0,   cos(rotate[1]) ], ]

        zmatrix = [ [  cos(rotate[2]),  sin(rotate[2]),  0 ],
                    [ -sin(rotate[2]),  cos(rotate[2]),  0 ],
                    [               0,               0,  1 ], ]

        ## frame * xmatrix * ymatrix * zmatrix
        ##
        ## But we'll do it in this order:
        ##  frame * ((xmatrix * ymatrix) * zmatrix )
        ##
        ##

        rmatrix = numarray.matrixmultiply( numarray.matrixmultiply( xmatrix, ymatrix ), zmatrix)
        frame = numarray.matrixmultiply(frame, rmatrix)

    ## Move the frame.
    if move != None:
        numarray.add(frame, move, frame)

    return frame
    

def tinkeratomreplace(infile, outfile, replacements):
    """Replace atom type numbers in a tinker xyz file.

    This function uses definitions in the second column (the string)
    to set atom type numbers (third column) in a tinker xyz file.
    It is designed to be used when Tinker's pdbxyz program can not
    select the right atom type numbers, for whatever reason.

    The first two arguments are the input .xyz file and the output  

    # all replacements must be a string of length 5, right justified !
    replacements = {}
    replacements["HT1 "] = "   88"
    replacements["HT2 "] = "   88"
    replacements["OT  "] = "  101"

    Returns the number of replacements made.

    {"OH2 ":"  186", "H1  ":"  187", "H2  ":"  187"}
    """
    oldxyz = file(infile)
    number_of_replacements = 0

    newxyz_list = []
    newxyz_list.append(oldxyz.readline())
    for line in oldxyz:
        atomname = line[8:12]
        if replacements.has_key(atomname):
            line = "".join([line[:48],  replacements[atomname],  line[53:]])
            newxyz_list.append( line )
            number_of_replacements += 1
        else:
            newxyz_list.append(line)
    newxyz_string = "".join(newxyz_list)
    #print newxyz
    outfo = file(outfile, "w")
    outfo.write(newxyz_string)
    outfo.close()
    print "number of replacements:", number_of_replacements
    return number_of_replacements


def _atomlist_to_leslist(atomlist, natoms):
    """Convert an "atomlist" to an "leslist"

    Originally used for locally-enhanced-sampling list creation.

    An atomlist lists each atom which is replicated.  An "leslist"
    lists the number of non-replicated atoms, the number of
    replicated atoms, the number of non-replicated atoms, etc.

    This function takes two arguments, "atomlist" is a list which
    contains an element specifying the number of each replicated,
    atom, and "natoms", the total number of atoms in the system.
    """
    type_ = False   # True = the last atom was in atomlist,
                    #        = the last atom was LESed
                    # False = the last atom wasn't in atomlist.
                    #        = the last atom was NOT LESed
    leslist = [ 0 ]
    for i in range(natoms):
        if i in atomlist:
            if type_ == True:
                leslist[-1] += 1
            else:
                type_ = True
                leslist.append(0)
                leslist[-1] += 1
        else:
            if type_ == False:
                leslist[-1] += 1
            else:
                type_ = False
                leslist.append(0)
                leslist[-1] += 1
    return leslist
                   
                
def normatomlist(atomlist):
    """Function to save time normalizing atom sets.

    Given an argument, this function will return it:
      - If the argument was a set, it will be converted to a list
        and sorted
      - If a list, it will be sorted
    """
    atomlist = list(atomlist)
    atomlist.sort()
    return atomlist

# Vector functions
_dot = numarray.dot
_norm = lambda x: math.sqrt(_dot(x,x))
_angle = lambda x, y: math.acos(_dot(x, y) / (_norm(x) * _norm(y)))


def _alignatomstoorigin(atoms):
    """

    atoms = ((a1, a1, a1),
             (a2, a2, a2),
             (a3, a3, a3))

    return (move, translist), the movement and rotations needed to bring
    about this transformation.
    """
    asarray = numarray.asarray
    M = asarray(atoms, type=numarray.Float32).copy()
    move = None      # What we are going to return
    translist = []   #

    # phase one: move to origin
    move = -M[0]
    M += move
    # rotate M1 to xz plane
    proj = ( M[1,0], M[1,1], 0 )
    ang = _angle(proj, (1,0,0))
    if proj[1] > 0:
        ang = -ang
    transform = (0, 0, ang)
    translist.append(transform)
    M = cordtransform(M, rotate=transform)
    
    print M[1,0], "should be zero- ahsuoe"
    # rotate M1 to along the x axis
    ang = _angle(M[1], (1,0,0))
    if M[1,2] < 0:
        ang = -ang
    transform = (0, ang, 0)
    translist.append(transform)
    M = cordtransform(M, rotate=transform)

    # bring last atom into x-y plane
    proj = (0, M[2,1], M[2,2])
    ang = _angle(proj, (0,1,0))
    if proj[2] < 0:
        ang = -ang
    transform = (-ang, 0, 0)
    translist.append(transform)
      #M = cordtransform(M, rotate=tranform)

    return (move, translist)


def alignframetoorigin(frame, atoms):
    """

    move frame to origin, aligning the atom indexes in "atoms" in the
    following manner:

    a1 is at the origin
    a2 is along the x axis
    a3 is in the xy plane
    """

    frame = frame.copy()
    move, translist = _alignatomstoorigin(atoms=(frame[atoms[0]],
                                                 frame[atoms[1]],
                                                 frame[atoms[2]]))
    frame += move

    for trans in translist:
        frame = cordtransform(frame, rotate=trans)
    return frame

def alignframetoposition(frame, atoms, position):
    """

    similar to alignframetoorigin, but aligns it to
    the place given by position=[vec1, vec2, vec3]
    """
    frame = frame.copy()
    frame = alignframetoorigin(frame, atoms)

    antimove, antitranslist = _alignatomstoorigin(atoms=position)
    antitranslist.reverse()
    antitranslist = [ -numarray.asarray(x) for x in antitranslist ]
    for antitrans in antitranslist:
        frame = cordtransform(frame, rotate=antitrans)
    frame -= antimove
    return frame


