
import mp3, sys, numarray


# Here is some argument parsing, but "optparse" (standard library in 2.3)
# is FAR better.
#
if len(sys.argv) < 3:    # do we have enough arguments?
        print "%s dcdname pdbname" % sys.argv[0]  # usage syntax
dcdname = sys.argv[1]   # get the file names
psfname = sys.argv[2]


# I always have trouble thinking of what the names will be.
# "mols" means "the molecules" we want to analyze
mols = mp3.system()               # set up the system here
mols.labels.getfrompsf(psfname)   # This opens the file, reads it in, closes it

mols.cord = mp3.dcd()             # mols.dcd isn't automatically there
                                  #  because there can be different types of "containers"
                                  #  for it
mols.cord.setfile(dcdname)        # Tells it where we are going to be reading from.
                                  #  You -no longer- need to call mols.cord.init()



group_of_interest = "OSPC"

print "Seaching for %s groups"%group_of_interest
print '---'

# find out which atoms we want to analyze
atoms_to_analyze = []      # a list to store them in
for i in xrange(0,mols.cord.natoms):   # loop over all atoms...
    if mols.labels.data[i].field('atomid') == group_of_interest:   # if it meets our arbitray criteria
        #print mols.labels.data[i]   # print out the line, (optional and slow)
        atoms_to_analyze.append(i)  # add that atom number to the list

# better:
atoms_to_analyze = [i for i in range(0,mols.cord.natoms) if mols.labels.data[i].field('atomid') == group_of_interest]

n_atoms_to_analyze = len(atoms_to_analyze)





#
# So, now, what analysis do we want to do ? Let's print out the average
#  x,y,z value for the first 5 frames.
#

for i in range(0,5):
    mols.cord.nextframe()
    xsum, ysum, zsum = 0, 0, 0   #tuple packing and unpacking !
    for atomnumber in atoms_to_analyze:
        xsum += mols.cord.frame[atomnumber,0]
        ysum += mols.cord.frame[atomnumber,1]
        zsum += mols.cord.frame[atomnumber,2]
    print "Frame %d: center: %f %f %f"%(i,xsum/n_atoms_to_analyze,ysum/n_atoms_to_analyze,zsum/n_atoms_to_analyze)

print "---"


#
# That way was rather slow. I'd prefer to use numarrays to do it all.
#
#
#
#

for i in range(5,10):
    mols.cord.nextframe()
    sum = numarray.zeros(shape=(3), type=numarray.Float32)       # "sum" is a numarray
    for atomnumber in atoms_to_analyze:  
        sum += mols.cord.frame[atomnumber]         #add a numarray to another numarray !
    print "Frame %d: center: %s"%(i, sum/n_atoms_to_analyze)
                                      #the numarray gets automatically converted to a string

#now wasn't THAT elegent ?
