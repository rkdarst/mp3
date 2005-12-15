#Richard Darst, 2005.
#
#

import numarray
import itertools
from math import pi
import mp3.log

from prdist_exclude import PDFExcluder

"""Pair Distribution utilities.

The key utilities for pair distribution functions:
- PDFHistogram, stores histograms.
- Cord objects output a histogram, to be added to other histograms.
- PDFExcluder, a class to aid in excluding certain interactions.

Keyword reference:
Below is a reference to the keywords needed for creating different
data structures relating to histograms.

Properties of histogram ():
* binwidth
* dim
- nbins (specify here, or implicitely when you first add bins)
- bins
- overflows
- volume

Properties needed to add points to a histogram:
* volume
* bins
* overflows

Properties needed to find a histogram frame:
* atomlist1
* atomlist2 (default to atomlist1)
* volume (???)
- boxsize (in order to wrap, required as of right now)
- nbins
"""

class PDFHistogram:
    def __init__(self, **keywords):
        """

        keywords (* = required, - = optional):
        
          * binwidth: the width of a bin.  Units are preserved
            throughout the calculation.

          * dim: How many dimensions does this PDF use?  should be an
            integer: 1, 2, or 3.

          - bins: the raw histogram bins.  The first bin should be [0,
            binwidth), the next [binwidth, binwidth*2),
            etc. (Optional, see 'addbins' method to incrementally add
            bins, or addpoints method to bin points.)  If you specify
            bins, you must specify overflows (probably zero)

          - overflows: see above.

          #- volume: Must be specified at least before bins are added.
          #  For a constant-total-volume system, specify it when the
          #  instance is created.  For a variable-volume boxsize,
          #  specify it when you do addbins().
          
          # ^ Depricated-- specify volume EVERY TIME now.  This makes
          # it more clear to follow.  You can specify it here if you
          # are specifing bins as well.
        
        """
        self._binwidth = keywords["binwidth"]
        self._dim = keywords["dim"]
        # volume and bins handled after we set up our dimensional consideration stuff.

        # function to return the amount of space in a bin.
        dr = self._binwidth
        dim = self._dim
        if dim == 1:
            # one dimensional -- it'l basically flat, but we'll do
            # this the "right" way anyway
            #
            # Area inside r == r
            # Area of slice == r+dr - r
            normvolfunc = lambda r: dr
        elif dim == 2:
            # Area inside r == pi*(r**2)
            # Area of slice == pi*((r+dr)**2) - pi*(r**2)
            normvolfunc = lambda r: pi*( ((r+dr)**2) - (r**2) )
        elif dim == 3:
            # Area inside r == (4/3)*pi*(r**3)
            # Area of slice == (4/3)*pi*((r+dr)**3) - (4/3)*pi*(r**3)
            normvolfunc = lambda r: (4./3)*pi*( (r+dr)**3 - (r**3) )
        elif dim >= 4:
            print "we can't handle more than three dimensions yet!"
            return None
        else:
            print "your dimension (%s) doesn't make sense"%dim
            return None
        self._normvolfunc = normvolfunc

        # we are given a constant volume for everything.
        #if keywords.has_key("volume"):
        #    self._volume = keywords["volume"]
        # error case: given bins but no volume.
        if keywords.has_key("bins") and not keywords.has_key("volume"):
            mp3.log.error("You must specify volumes before you specify bins.")
        # we are given bins to initilize with.
        if keywords.has_key("bins"):
            self.addbins(bins=keywords["bins"],
                         overflows=keywords["overflows"],
                         volume=keywords["volume"])
        # We don't want to add any bins to start off with, but we want
        # to get it set up for adding just points later on.
        if keywords.has_key("nbins"):
            nbins = keywords["nbins"]
            self._bins = numarray.zeros(nbins)
            self._volbins = numarray.zeros(nbins, type=numarray.Float32)
            self._overflows = 0


    def __iadd__(self, other):
        """Add two PDFHistograms together, checking that key attributes match.
        
        This should be used to do in-place addition of two Pdfhistograms (since
        we're basically doing in-place modification when adding bins, anyway.)
        For example, each Cord.pairdistfunc() returns a PDFHistogram.  So when 
        you are accumuliting a pair distribution function over many frames, you
        could do something like:

        H += Cord.pairdistfunc()
        """
        if len(self._bins) != len(other._bins):
            mp3.log.error("You must have the number of bins match in both histograms!")
            raise "" #FIXME
        if self._binwidth != other._binwidth:
            mp3.log.error("You must have equal binwidths both histograms!")
            raise "" #FIXME

        #xor = lambda a,b: (a or b) and not (a and b)
        #if xor(hasattr(self, "_volume"), hasattr(other, "_volume")):
        #    mp3.log.error("Both histograms must have a constant volume or not!")
        #    #raise "" #FIXME

        self._overflows += other._overflows
        self._bins += other._bins
        self._volbins += other._volbins
        
        
        # if hasattr(self, "_volume"):
        #     H = PDFHistogram()
        return self

    def __isub__(self, other):
        """Subtract in place.

        Exactly the opposite of add-in-place.
        """
        if len(self._bins) != len(other._bins):
            mp3.log.error("You must have the number of bins match in both histograms!")
            raise "" #FIXME
        if self._binwidth != other._binwidth:
            mp3.log.error("You must have equal binwidths both histograms!")
            raise "" #FIXME

        self._overflows -= other._overflows
        self._bins -= other._bins
        self._volbins -= other._volbins

        return self

    #
    # Putting information into the object.
    #

    def addbins(self, bins, overflows, volume ):
        """Add more data to this histogram.

        Add bins incrementally.  The argument 'bins' must be an array
        of binned data, with the proper length.  'overflows' should be
        how many points are too far away to be binned.  The 'volume'
        keyword should be the volume of this frame, if it hasn't
        already been specified.
        """
        # In general, we have a problem here.  We can add data to the
        # histogram at any time, but we also have post-processing that
        # we have done to it!  This other data must be updated.  We
        # can either generate this other data dynamically, or have
        # some way of updating it.

        normvolfunc = self._normvolfunc
        binwidth = self._binwidth

        #if volume is None:
        #    # if volume not specified, it must already be known.
        #    volume = self._volume
            
        # convert inputs into numarrays
        bins = numarray.asarray(bins)
        volbins = numarray.asarray(bins, type=numarray.Float32 )
        # do volume corrections
        for bin, i in zip(volbins, itertools.count()):
            r = binwidth * i
            correction = volume / ( normvolfunc(r))
            volbins[i] = bin * correction
        volbins[0] = 0

        # store data.
        if not hasattr(self, "_bins"):
            self._bins = bins
            self._volbins = volbins
            self._overflows = overflows
        else:
            self._bins = numarray.add(self._bins, bins)
            self._volbins = numarray.add(self._volbins, volbins)
            self._overflows += overflows

    def addpoints(self, points, volume):
        """Add points from a distance list.

        This function will bin points for you add them using the
        add_bins method.  It knows the binwidth based on the binwidth
        used to initialize the histogram class.

        However, if we haven't added any points yet, it won't know how
        many bins there are.  In this case, we must initilize the
        PDFHistogram with a 'nbins=' paramater, so that at this point
        we'll know how far out to go.

        """
        # get the information needed to do binning.
        nbins = len(self._bins)
        binwidth = self._binwidth
        max_distance = nbins*binwidth
        # set up variables.
        overflows = 0
        bins = [0] * nbins
        # do binning
        for distance in points:
            if distance >= max_distance:
                overflows += 1
            else:
                bins[int(distance/binwidth)] += 1
        # add in results.
        self.addbins(bins, overflows, volume)

    

    #
    # Getting information out of the object.
    # 

    def total_count(self):
        """Count of total number of samplings.

        Return total number of points in the bins.

        NOTE: this includes zero-distance overlaps, when atomN was
        correllated with atomN (zero-distance pairs.  See the
        'total_pdf_count' method for a count not including these
        overlaps.
        """
        self._total_bin_count = self._bins.sum() + self._overflows
        return self._total_bin_count
    def total_pdf_count(self):
        """Count of non-overlaping samplings.
        """
        self._total_pdf_count = self._bins[1:].sum() + self._overflows
        return self._total_pdf_count

    def bins_d(self):
        """Distances of each bin.

        Returns the lower bound of each bin."""
        return self._binwidth * numarray.arange(len(self._bins))

    def bins(self):
        """ Return unnormalized bins.
        """
        return self._bins

    def normalized_bins(self):
        """Return (bins / volume) / total_pdf_count
        """
        total_pdf_count = self.total_pdf_count()
        return numarray.asarray(self._volbins / total_pdf_count)

    def integrated_bins(self):
        n_bins = self.normalized_bins()
        i_bins = [ ]
        for i in range(len(n_bins)):
            i_bins.append(sum(n_bins[0:i])*self._binwidth )
        return i_bins

def pdftrj(**keywords):
    """Take a PDF of a whole cord object.

    This function will take a cord object (trajectory of frames) and
    find the PDFHistogram of each frame in it.

    Most keywords here are passed through unchanged to the
    PDFHistogram class, but we have some extra considerations:

    keyword arguments:

      * Any keywords needed for the Cord.pairdistfunc method
        (atomlist1, [atomlist2], binwidth, nbins, dim)

      * cord: the cord object to take the frames of.  We call
        nextframe() _before_ we take each histogram, so it should be
        set to the frame before the one we want to begin at.  Note
        that if you want to do all frames, the cord object begins an
        the right place.

      * boxsize: If it's a boxsize object (if it's callable), it will
        be used to read off successize boxsizes.  If it's not
        callable, it's a tuple-like object used for the constant box
        size.

      - nframes: Optional.  If it's given, 
    """
    cord = keywords["cord"]
    if keywords.has_key("nframes"):
        nframes = keywords["nframes"]
    else:
        nframes = cord.nframes()

    #cord.nextframe()
    H = cord.pairdistfunc(**keywords)  # use bins to initilize it, using nbins.
    for i in range(nframes):
        cord.nextframe()
        H += cord.pairdistfunc(**keywords)

    return H

            

def _plot(bins, bins_d):
    import Gnuplot
    g = Gnuplot.Gnuplot()
#    g("set yrange [:3]")
    array = zip(bins_d, bins)
    item = Gnuplot.Data(array, with="lines")
    g.plot(item)
    raw_input()

def cpdf_histframe(**keywords):
    """Thin wrapper around the C PDF histrogram function.
    """
    keywords["frame"] = numarray.asarray(keywords["frame"])
    keywords["atomlist1"] = tuple(keywords["atomlist1"])
    if not hasattr(keywords, "atomlist2"):
        keywords["atomlist2"] = keywords["atomlist1"]
    else:
        keywords["atomlist2"] = tuple(keywords["atomlist2"])
    keywords["binwidth"] = float(keywords["binwidth"])
    keywords["boxsize"] = tuple( [ float(keywords["boxsize"][i]) for i in range(3) ] ) 
    keywords["nbins"] = int(keywords["nbins"])

    import sys
    sys.path.append("/home/richard/research/mp3/src/mp3c")
    import cpdf
    return_dict = cpdf.cpdf_histframe( (keywords,) )
    #return None
    # return value has bins, overflows. in it.
    return return_dict
