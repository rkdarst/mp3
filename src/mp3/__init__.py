#  MP3 -- Molecular Processing Package
#
#  This package aims to provide tools for manipulating common molecular
# file formats, and providing a method of analysis.
#
#  Developed at the University of Texas at Austin
#  Contact: (rkd zgib.net)
#


import logging, sys
_log = logging.getLogger('mp3')
_theloghandler = logging.StreamHandler(sys.stdout)
_theloghandler.setFormatter( logging.Formatter("mp3:%(levelname)s:%(filename)s %(message)s"))
_log.addHandler(_theloghandler)
_log.debug('Logger Initilized')
_log.info('DCD Management Library 0.0')
_log.debug('Old joke: never use version 0.x or x.0 (UT ITS UNIX Group)')
#thelog.setLevel(logging.INFO)

#  The purpose of the __init__ module is to load my definitions from all
# other files.  I hope that this will make it more organized.
#

from labels import *         # Labels manages, such as atomnames and types
from system import *         # Combination of labels, bonds, coordinates, ...
from cord import *           # Parent class for coordinate objects
from cordalign import *      # Aligns atoms in a cord
from cordatomslice import *  # Extract only some atoms from cord object
from cordcenterer import *   # Translate coordinates to recenter system
from corddcd import *        # Reads DCD files
from corddummy import *      # Example of making a wrapper for cords
from cordlesmod import *     # Wrapper for locally-enhanced-sampling cords.
from cordmerge import *      # Merges multiple coordinates into a continous set
from cordminimage import *   # Wraps coordinates to their minimum images
from cordpdb import *        # Reads coordinates from PDBs
from cordtinkerxyz import *  # Reads coordinates from Tinker XYZ files.
from cordtinkerarc import *  # Reads coordinates from Tinker archive files.
from cordtransform import *  # Cord transformations (move+rotate)
from functions import *      # Useful functions
from msder import *          # A class to aid with mean-square-deviation analyzes.
from xst import *            # Deals with NAMD xst files (periodic box size)


#from xyz1 import *          # Incomplete module to read xyz coordinates
