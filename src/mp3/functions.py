# This module contains functions to simplify various processing
# using mp3.
#
#
#

import numarray, math
import mp3

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
