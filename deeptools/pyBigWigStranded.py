#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re

import numpy as np
import pyBigWig

class pyBigWig1:
    """Wrapper around pyBigWig to ignore strand.

    This modifies some pyBigWig methods to expect (and ignore) strand
    information. This way, callers can get data either from
    a pyBigWig object, or a strandedPyBigWig object, without
    knowing whether the file is stranded.

    Possibly this module belongs in pyBigWig.

    Also, this module (currently) requires numpy support in pyBigWig.
    """
    def __init__(self, bigWigFile):
        """Constructor is a wrapper around pyBigWig.pyBigWig.
        """
        self.bigWig = pyBigWig.open(bigWigFile)

    def values(self, chrom, start, end, strand='+'):
        """Gets values at a region.

        This is like pyBigWig.values(), except that expects the strand argument
          (which it ignores).
        ??? possibly this should reverse values, if strand is '-' ?
        """
        return self.bigWig.values(chrom, start, end)

    def chroms(self, chroms=None):
        """This is just forwarded.
        """
        if chroms:
            return self.bigWig.chroms(chroms)
        else:
            return self.bigWig.chroms()


class strandedPyBigWig:
    """Wrapper for accessing pairs of 'stranded' bigWig files.

    This has almost the same API as pyBigWig, except that
    most methods include a 'strand' argument. When these are
    called, this gets the result from the appropriate stranded file.

    Note that that these methods return the absolute value of the
    contents of the bigWig. This is because the - strand file is
    sometimes stored negated. It seems simpler to just always
    return the absolute value, than to add a flag indicating
    whether the - strand file is stored negated.
    """

    def __init__(self, plusStrandFile,
        suffixPairs=['fwd.bw,rev.bw', 'plus.bw,minus.bw', 'pl.bw,mn.bw', 'p.bw,m.bw'],
        minusStrandFile=None):
        """Constructor opens a pair of bigWig files for reading.

        suffixPairs: these define the + and - filename pairing.
        minusStrandFile: the minus strand file can also simply be given
          explicitly (this overrides the name computed using suffixPairs)
        """
        # if - strand file isn't given, attempt to compute it
        if not minusStrandFile:
            minusStrandFile = getMinusStrandFile(plusStrandFile, suffixPairs)
        # if we still don't have a - strand filename, throw an error
        if not minusStrandFile:
            raise ValueError('couldn\'t compute - strand filename')
        # open files
        self.plusBigWig = pyBigWig1(plusStrandFile)
        self.minusBigWig = pyBigWig1(minusStrandFile)

    def chroms(self, chroms=None):
        """Gets chromosomes.

        This just uses the + strand file's chromosomes (without
        checking that they agree with the - strand file's).
        """
        if chroms:
            return self.plusBigWig.chroms(chroms)
        else:
            return self.plusBigWig.chroms()

    def values(self, chrom, start, end, strand='+'):
        """Gets values at a region specified by an object.

        This will get values from plusBigWig or minusBigWig,
          depending on the strand.
        Returns: absolute value of values at that region.
          If strand is '-', these will be in descending coordinate
          order; otherwise they'll be in ascending coordinate order.
        """
        # get the relevant bigWig object
        bigWig = self.minusBigWig if strand == '-' else self.plusBigWig
        # forward this call to that object, taking absolute value
        return np.abs(bigWig.values(chrom, start, end, strand))

def getMinusStrandFile(plusStrandFile,
        suffixPairs=['fwd.bw,rev.bw', 'plus.bw,minus.bw', 'pl.bw,mn.bw', 'p.bw,m.bw']):
    """Converts the + strand filename to the - strand filename.

    plusStrandFilename: name of the + strand file
    suffixPairs: these define the + and - filename pairing.
      For each A,B pair, if plusStrandFile ends with A, then minusStrandFile is
      assumed to end with B.
    Returns: name of the - strand file, or None if the + strand file's
      name didn't match any of the known + strand file's suffixes.
    """
    # loop through + suffix, until a matching one is found
    for suffixPair in suffixPairs:
        (plusSuffix, minusSuffix) = suffixPair.split(',')
        # try replacing the suffix (at the end of the filename)
        minusStrandFile = re.sub(re.escape(plusSuffix) + '$', minusSuffix,
            plusStrandFile)
        # if this is different, then the replacement happened
        if minusStrandFile != plusStrandFile:
            return minusStrandFile
    # at this point, none of the suffixes matched
    return None

def openBigWigGuessingStrandedness(bigWigFile,
        suffixPairs=['fwd.bw,rev.bw', 'plus.bw,minus.bw', 'pl.bw,mn.bw', 'p.bw,m.bw']):
    """Opens a bigWig file/files, guessing whether they're stranded.

    bigWigFile: name of the file (or the + strand file)
    suffixPairs: these define the  and - filename pairing
      (as elsehwere)
    Returns: an object for reading from a bigWig file, either stranded
      (if bigWigFile ends with one of the suffixes), or not
      (if it doesn't).
    """
    # see if bigWigFile "looks" stranded
    minusStrandFile = getMinusStrandFile(bigWigFile, suffixPairs)
    if minusStrandFile:
        # if so, assume it's a stranded pair
        return strandedPyBigWig(bigWigFile, minusStrandFile=minusStrandFile)
    else:
        # otherwise, assume it's unstranded
        return pyBigWig1(bigWigFile, suffixPairs)
