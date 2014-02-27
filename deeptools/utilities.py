import sys
import os

debug = 0


def getGC_content(dnaString, as_fraction=True):
    if len(dnaString) == 0:
        return None
    if dnaString.count('N') > len(dnaString) * 0.05:
        raise Exception("too many NNNs in assembly sequence")
        return None

    gc = 0
    gc += dnaString.count('G')
    gc += dnaString.count('g')
    gc += dnaString.count('C')
    gc += dnaString.count('c')
    if as_fraction:
        return(float(gc) /len(dnaString))
    else:
        return gc


def tbitToBamChrName(tbitNames, bamNames):
    """ checks if the chromosome names from the two-bit and bam file coincide.
        In case they do not coincide, a fix is tried. If successful, then
        a mapping table is returned.
        tbitNames and bamNames should be lists
    """

    chrNameBitToBam = dict((x, x) for x in tbitNames)
    if set(bamNames) != set(tbitNames):
        sys.stderr.write("Bam and 2bit do not have matching "
                         "chromosome names:\n2bit:{}\n\nbam:{}"
                         "\n\n".format(tbitNames, bamNames))

        if set(["chr" + x if x != 'dmel_mitochondrion_genome'
                else 'chrM' for x in bamNames]) == set(tbitNames):
            sys.stderr.write("Adding 'chr' seems to solve the problem. "
                             "Continuing ...")
            chrNameBitToBam = dict([("chr" + x
                                     if x != 'dmel_mitochondrion_genome'
                                     else 'chrM', x) for x in bamNames ])
        elif set([x for x in tbitNames if x.count('random') == 0
                  and x.count('chrM') == 0]) == set(bamNames):
            if debug:
                print "Removing random and mitochondrial chromosomes"\
                    "fixes the problem"
            chrNameBitToBam = dict([(x, x) for x in tbitNames
                                    if x.count('random') == 0 and
                                    x.count('chrM') == 0])
        elif len(set(bamNames).intersection(set(tbitNames)) ) > 0:
            if debug:
                print "Using only common chromosomes between between "
                "index and reference:"
                print set(bamNames).intersection(set(tbitNames))
            chrNameBitToBam = dict([(x, x) for x in
                                    set(bamNames).intersection(set(tbitNames))])
        else:
            if debug:
                print "Index and reference do not have matching "
                "chromosome names"
            exit(0)

    return chrNameBitToBam


def getCommonChrNames(bamFileHandlers, verbose=True):
    r"""
    Compares the names and lengths of a list of bam file handlers.
    The input is list of pysam file handlers.

    The function returns a duple containing the common chromosome names
    and the common chromome lengths.

    Hopefully, only _random and chrM are not common.
    """
    outMessage = []
    commonChr = set( [ "{}|{}".format(bamFileHandlers[0].references[i], bamFileHandlers[0].lengths[i]) for i in range(0,len(bamFileHandlers[0].references)) ] )
    
    maxRefLen = len(commonChr)
    for j in range(1, len(bamFileHandlers)):
       refLen = len(bamFileHandlers[j].references)
       commonChr = commonChr &  set ([ "{}|{}".format(bamFileHandlers[j].references[i], bamFileHandlers[j].lengths[i]) for i in range(0,refLen) ] )
       if refLen > maxRefLen:
          maxRefLen = refLen
       
    if len(commonChr) != maxRefLen:
        outMessage.append( "\nReferenced chromosome names in the given bam files differ.\n\n" )
        for i in range(0, len(bamFileHandlers)):
           outMessage.append( "{0:>15}\t{1:>10}".format("chr names bam {}".format(i+1),"chr length") )
        outMessage.append( "\n\n" )

        # get the largest number of references (i.e. chromosome names)
        # from the list of bam files
        # maxI = max( map( (lambda x: len(x.references)), bamFileHandlers ) )
        maxI = max( [ len(x.references) for x in bamFileHandlers ] )
        for i in range(0, maxI):
           for j in range(0, len(bamFileHandlers)):
              try:
                 outMessage.append( "{0:>15}\t{1:>10}".format(bamFileHandlers[j].references[i], 
                                                 bamFileHandlers[j].lengths[i]) )
              except:
                 outMessage.append( "{0:^15}\t{1:^10}\t".format("--", "--") )
           outMessage.append( "\n" ) # force a new line

    if len(commonChr) == 0:
       outMessage.append( "No common chromosomes found." )
       outMessage.append( "\nAre these bam files from different species or different assemblies?\n" )
       sys.stderr.write("".join(outMessage) )
       exit(1)

    # the common chromosomes has to be sorted as in the original
    # bam files
    chrSizes = []
    for i in range(0, len(bamFileHandlers[0].references)):
        if "{}|{}".format(bamFileHandlers[0].references[i],
                          bamFileHandlers[0].lengths[i]) in commonChr:

            chrSizes.append((bamFileHandlers[0].references[i],
                             bamFileHandlers[0].lengths[i]))

    outMessage.append("\nUsing the following set of common chromosome "
                      "names and lengths:\n")
    for chrSize in chrSizes:
        outMessage.append("{0:>15}\t{1:>10}\n".format(chrSize[0], chrSize[1]))

    if verbose:
        sys.stderr.write("".join(outMessage))
    return chrSizes


def copyFileInMemory(filePath, suffix=''):
    """
    copies a file into the special /dev/shm device which
    moves the file into memory.
    This process speeds ups the multiprocessor access to such files
    """

    # fallback for windows users
    if os.name == 'nt':
        return filePath

    memFileName = getTempFileName(suffix=suffix)
    import shutil
    shutil.copyfile(filePath, memFileName)

    return memFileName


def getTempFileName(suffix=''):
    """
    returns a temporary file name.
    If the special /dev/shm device is available,
    the temporary file would be located in that folder.
    /dv/shm is a folder that resides in memory and
    which has much faster accession.
    """
    import tempfile
    # is /dev/shm available?
    try:
        _tempFile = tempfile.NamedTemporaryFile(prefix="_deeptools_",
                                                suffix=suffix,
                                                dir='/dev/shm',
                                                delete=False)
    except OSError:
        _tempFile = tempfile.NamedTemporaryFile(suffix=suffix,
                                                delete=False)

    memFileName = _tempFile.name
    _tempFile.close()
    return memFileName


def which(program):
    """ method to identify if a program
    is on the user PATH variable.
    From: http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
