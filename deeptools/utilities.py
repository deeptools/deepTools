import sys
import os

debug = 0


def getGC_content(dnaString, as_fraction=True):
    if len(dnaString) == 0:
        return None
    if dnaString.count('N') > len(dnaString) * 0.05:
        raise Exception("WARNING: too many NNNs present in sequence of length {}".format(len(dnaString)))
        return None

    gc = 0
    gc += dnaString.count('G')
    gc += dnaString.count('g')
    gc += dnaString.count('C')
    gc += dnaString.count('c')
    if as_fraction:
        return(float(gc) / len(dnaString))
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

        if len(set(bamNames).intersection(set(tbitNames))) > 0:
            sys.stderr.write("Using the following common chromosomes between "
                             "bam chromosome names and 2bit chromosome "
                             "names:\n")
            for item in set(bamNames).intersection(set(tbitNames)):
                sys.stderr.write(item + "\n")
            chrNameBitToBam = dict([(x, x) for x in
                                    set(bamNames).intersection(set(tbitNames))])
        elif set(["chr" + x if x != 'dmel_mitochondrion_genome'
                  else 'chrM' for x in bamNames]) == set(tbitNames):
            sys.stderr.write("Adding 'chr' seems to solve the problem. "
                             "Continuing ...")
            chrNameBitToBam = dict([("chr" + x
                                     if x != 'dmel_mitochondrion_genome'
                                     else 'chrM', x) for x in bamNames])
        elif set([x for x in tbitNames if x.count('random') == 0 and
                 x.count('chrM') == 0]) == set(bamNames):
            if debug:
                print "Removing random and mitochondrial chromosomes"\
                    "fixes the problem"
            chrNameBitToBam = dict([(x, x) for x in tbitNames
                                    if x.count('random') == 0 and
                                    x.count('chrM') == 0])
        elif len(set(["chr" + x for x in bamNames if x != 'dmel_mitochondrion_genome']).intersection(set(tbitNames))) > 0:
            bamNames2 = ["chr" + x for x in bamNames if x != 'dmel_mitochondrion_genome']
            sys.stderr.write("Adding 'chr' seems to solve the problem for the following "
                             "chromosomes...")
            for item in set(bamNames2).intersection(set(tbitNames)):
                sys.stderr.write(item + "\n")

            chrNameBitToBam = {"chrM": "MT"}
            for i in xrange(len(bamNames)):
                if bamNames2[i] in tbitNames:
                    chrNameBitToBam.update({bamNames2[i]: bamNames[i]})
        elif len(set([x[3:] for x in bamNames if x.startswith("chr")]).intersection(set(tbitNames))) > 0:
            bamNames = [x for x in bamNames]
            bamNames2 = [x[3:] for x in bamNames if x.startswith("chr")]
            if debug:
                sys.stderr.write("Removing 'chr' seems to solve the problem for the following "
                                 "chromosomes...")
                for item in set(bamNames).intersection(set(tbitNames)):
                    sys.stderr.write(item + "\n")

            chrNameBitToBam = {"MT": "chrM"}
            for i in xrange(len(bamNames)):
                if bamNames2[i] in tbitNames:
                    chrNameBitToBam.update({bamNames2[i]: bamNames[i]})
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
    def get_chrom_and_size(bam_handler):
        """
        Reads the chromosome/scaffold name and the length from
        the bam file and returns a list of (chromname, size) tuples
        :param bam_handler:
        :return: list of (chrom, size) tuples
        """
        return [(bam_handler.references[i], bam_handler.lengths[i])
                for i in range(0, len(bam_handler.references))]

    def print_chr_names_and_size(chr_set):
        sys.stderr.write("chromosome\tlength\n")
        for name, size in chr_set:
            sys.stderr.write("{0:>15}\t{1:>10}\n".format(name, size))

    common_chr = set(get_chrom_and_size(bamFileHandlers[0]))
    non_common_chr = set()

    for j in range(1, len(bamFileHandlers)):
        _names_and_size = set(get_chrom_and_size(bamFileHandlers[j]))
        if len(common_chr & _names_and_size) == 0:
            #  try to add remove 'chr' from the chromosome name
            _corr_names_size = set()
            for chrom_name, size in _names_and_size:
                if chrom_name.startswith('chr'):
                    _corr_names_size.add((chrom_name[3:], size))
                else:
                    _corr_names_size.add(('chr' + chrom_name, size))
            if len(common_chr & _corr_names_size) == 0:
                message = "No common chromosomes found. Are the bam files files " \
                          "from the same species and same assemblies?\n"
                sys.stderr.write(message)
                print_chr_names_and_size(common_chr)

                sys.stderr.write("\nand the following is the list of the unmatched chromosome and chromosome\n"
                                 "lengths from file\n{}\n".format(bamFileHandlers.name))
                print_chr_names_and_size(_names_and_size)
                exit(1)
            else:
                _names_and_size = _corr_names_size

        non_common_chr |= common_chr ^ _names_and_size
        common_chr = common_chr & _names_and_size

    if len(non_common_chr) > 0:
        sys.stderr.write("\nThe following chromosome names did not match between the the bam files\n")
        print_chr_names_and_size(non_common_chr)

    # the common chromosomes has to be sorted as in the original
    # bam files
    chr_sizes = []
    for tuple in get_chrom_and_size(bamFileHandlers[0]):
        if tuple in common_chr:
            chr_sizes.append(tuple)

    return chr_sizes, non_common_chr


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
    import config as cfg
    # get temp dir from configuration file
    tmp_dir = cfg.config.get('general', 'tmp_dir')
    if tmp_dir == 'default':
        _tempFile = tempfile.NamedTemporaryFile(prefix="_deeptools_",
                                                suffix=suffix,
                                                delete=False)

    else:
        try:
            _tempFile = tempfile.NamedTemporaryFile(prefix="_deeptools_",
                                                    suffix=suffix,
                                                    dir=tmp_dir,
                                                    delete=False)
        # fall back to system tmp file
        except OSError:
            _tempFile = tempfile.NamedTemporaryFile(prefix="_deeptools_",
                                                    suffix=suffix,
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
