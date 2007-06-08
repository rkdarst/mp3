# Richard Darst, 2005

import numpy

import mp3

class PDFExcluder:
    """Produce histograms of complicated atomlists.

    This class is used to exclude certain interactions when finding
    pair distribution functions.  There is a C module for finding
    histograms, but it does not handle exclusion of certain
    interactions.  In order to do that, the present method is to find
    the histogram of all interactions, find the histogram of only the
    interactions which are to be excluded, and subtract them.

    #FIXME undocumented.

    Initilization arguments:
    * atomlist1
    * atomlist2
    * bonds : the list of interactions to consider
    * **keywords - keywords needed when creating the histogram for each frame.

    assumptions: don't reorder atomlists.

    """
    def __init__(self, atomlist1, atomlist2, bonds, **keywords):
        atomlist1 = atomlist1
        atomlist2 = atomlist2
        self._keywords = keywords
        # the primary loop is done over atomlist1 We want this to be
        # the short loop.

        # NOTE: don't let it the lists reordered

        # bondedexcludelist[i] is which bonds to subtract from
        # atomlist1[i]

        excludelist = [ ]
        for atom in atomlist1:
            excludelist.append( [i for i in bonds[atom] if i in atomlist2 ] )
        print excludelist[:10]

        self._atomlist1 = atomlist1
        self._excludelist = excludelist

    def excludedhistogram(self, cord, volume):
        """Return a histogram of excluded bins

        Arguments:

        * cord: the cord object to process
        * volume: the volume, same definition as for histograms.

        #FIXME undocumented
        """
        #
        # Notice
        #
        
        # This code was taken from cord.py, Cord.pairdistfunc.  I have
        # duplicated it here, instead of making it modular (yet) for
        # efficiency reasons.  Both parts should be updated in unison.

        frame = cord.frame()
        keywords = self._keywords
        excludelist = self._excludelist

        binwidth = keywords["binwidth"]
        nbins = keywords["nbins"]
        boxsize = keywords["boxsize"]
        dim = keywords["dim"]
        # set atomlists, we default to using the same atomlist for both.

        # find volume
        #if not keywords.has_key("volume"):
        #    volume = boxsize[0]*boxsize[1]*boxsize[2]
        #else:
        #    volume = keywords["volume"]

        # do we want to use C?
        useC = False

        # This loop choses if we should do the binning in C or Python.
        if useC:
            pass  # this has been excluded above
            ## Use the C module to do the binning.
            #keywords["frame"] = frame
            #keywords["atomlist1"] = atomlist1
            #keywords["atomlist2"] = atomlist2
            #return_dict = mp3.prdistfunc.cpdf_histframe(**keywords)
            ##return None
            #bins = return_dict["bins"]
            #overflows = return_dict["overflows"]

        else:
            # Do binning in Python.
            bins = [0] * nbins
            overflows = 0                # how many points are too large to count.
            max_distance = nbins*binwidth # distance => this counts as an overflow
            
            # we use numpy array to loop over atomlist 2.  This saves us a lot of time.
            for i, atom1 in enumerate(self._atomlist1):
                atomlist2 = excludelist[i]
                deltas = frame[atomlist2] - frame[atom1]   # displacements
                numpy.divide(deltas, boxsize, deltas)   # new coordinates:
                                                           # [0,boxsize)  ->  [0,1)
                # These lines map all points to the interval [-.5, .5)
                shift = numpy.add(deltas, .5)           
                numpy.floor(shift, shift)
                numpy.subtract(deltas, shift , deltas)
                numpy.multiply(deltas, boxsize, deltas)  # --> [-size/2, size/2 )
                deltas *= deltas
                deltas = numpy.sum(deltas, axis=1)
                #print deltas
                distances = numpy.sqrt(deltas)
                
                # record the data
                for distance in distances:
                    if distance >= max_distance:
                        overflows += 1
                    else:

                        bins[int(distance/binwidth)] += 1


        # build the return histogram
        H = mp3.prdistfunc.PDFHistogram(bins=bins,
                                        overflows=overflows,
                                        volume=volume,
                                        binwidth=binwidth,
                                        dim=dim)
        #H.addbins(bins, overflows=overflows)
        return H
                
