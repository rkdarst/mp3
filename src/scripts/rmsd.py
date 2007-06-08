
import mp3, sys, numpy
from math import sqrt


if len(sys.argv)< 3:    # do we have enough arguments?
    print "%s dcdname pdbname" % sys.argv[0]  # usage syntax
dcdname = sys.argv[1]   # get the file names
pdbname = sys.argv[2]


mols = mp3.system()    # set up the system here
mols.labels.getfrompsf(pdbname)
mols.cord = mp3.dcd()
mols.cord.setfile(dcdname)   #this will init it now


# find out which atoms we want to analyze
atoms_to_analyze = []      # a list to store them in
for i in xrange(0,mols.cord.natoms):   # loop over all atoms...
    if mols.labels.data[i].field('atomid') == 'OSPC':   # if it meets our arbitray criteria
        print mols.labels.data[i]   # print out the line, (optional and slow)
        atoms_to_analyze.append(i)  # add that atom number to the list


#print atoms_to_analyze   # more verification of what we are going to analyze
natoms_to_analyze = len(atoms_to_analyze)   # hom many atoms are we analyzing here ?

c_win = 5   # how long we want the MSD window to be
data = numpy.zeros(dtype=numpy.float32, shape=(natoms_to_analyze, c_win, 3))
    #this array will hold the last c_win frames
    # data[n,t,3] , where n is the atom, and t how new this data is. current frame, t=0
    
msd = numpy.zeros(dtype=numpy.float32, shape=(natoms_to_analyze, c_win))
    #this array holds the sum of the squares of the displacements
    
ntrials = 0  # we need to keep track of how many samples we get over time

for i in xrange(c_win-1, -1, -1):   # go backwards
    mols.cord.nextframe()          # getting the next frame
    mols.cord.frame[i,j]
    data[:,i] = mols.cord.frame[atoms_to_analyze]  #and storing it in data
                                                 #notice how we slice it with a -list- !


for ii in xrange(c_win,mols.cord.nframes):         # until we are through with this dcd.
    for i in xrange(0,c_win):           # do a test with all of the date we have stored
        msd[:,i] += ( numpy.sum ((data[:,0,:]-data[:,i,:]) * (data[:,0,:]-data[:,i,:]), axis=1))
    ntrials += 1
    data[1:] = data[:-1]
    data[:,-1,:] = mols.cord.nextframe()[atoms_to_analyze]


numpy.divide( msd, ntrials, msd)
print msd[0:10]

