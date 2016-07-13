import sys
import os
import pysam
from deeptoolsintervals import GTF
from deeptools.bamHandler import openBam


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
                print("Removing random and mitochondrial chromosomes"
                      "fixes the problem")
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
            for i in range(len(bamNames)):
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
            for i in range(len(bamNames)):
                if bamNames2[i] in tbitNames:
                    chrNameBitToBam.update({bamNames2[i]: bamNames[i]})
        else:
            if debug:
                print("Index and reference do not have matching ")
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
    from deeptools import config as cfg
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


def gtfOptions(allArgs=None):
    """
    This is used a couple places to setup arguments to mapReduce
    """
    transcriptID = "transcript"
    exonID = "exon"
    transcript_id_designator = "transcript_id"
    keepExons = False
    if allArgs is not None:
        allArgs = vars(allArgs)
        transcriptID = allArgs.get("transcriptID", transcriptID)
        exonID = allArgs.get("exonID", exonID)
        transcript_id_designator = allArgs.get("transcript_id_designator", transcript_id_designator)
        keepExons = allArgs.get("keepExons", keepExons)
    return transcriptID, exonID, transcript_id_designator, keepExons


def toString(s):
    """
    This takes care of python2/3 differences
    """
    if isinstance(s, str):
        return s
    if isinstance(s, bytes):
        if sys.version_info[0] == 2:
            return str(s)
        return s.decode('ascii')
    if isinstance(s, list):
        return [toString(x) for x in s]
    return s


def toBytes(s):
    """
    Like toString, but for functions requiring bytes in python3
    """
    if sys.version_info[0] == 2:
        return s
    if isinstance(s, bytes):
        return s
    if isinstance(s, str):
        return bytes(s, 'ascii')
    if isinstance(s, list):
        return [toBytes(x) for x in s]
    return s


def mungeChromosome(chrom, chromList):
    """
    A generic chromosome munging function. "chrom" is munged by adding/removing "chr" such that it appears in chromList

    On error, None is returned, but a common chromosome list should be used beforehand to avoid this possibility
    """
    if chrom in chromList:
        return chrom

    if chrom == "MT" and "chrM" in chromList:
        return "chrM"
    if chrom == "chrM" and "MT" in chromList:
        return "MT"

    if chrom.startswith("chr") and chrom[3:] in chromList:
        return chrom[3:]
    if "chr" + chrom in chromList:
        return "chr" + chrom

    # This shouldn't actually happen
    return None


def bam_total_reads(bam_handle, chroms_to_ignore):
    """Count the total number of mapped reads in a BAM file, filtering
    the chromosome given in chroms_to_ignore list
    """
    if chroms_to_ignore:
        import pysam

        lines = pysam.idxstats(bam_handle.filename)
        lines = toString(lines)
        if type(lines) is str:
            lines = lines.strip().split('\n')
        if len(lines) == 0:
            # check if this is a test running under nose
            # in which case it will fail.
            if len([val for val in sys.modules.keys() if val.find("nose") >= 0]):
                sys.stderr.write("To run this code inside a test use disable "
                                 "output buffering `nosetest -s`\n".format(bam_handle.filename))
            else:
                sys.stderr.write("Error running idxstats on {}\n".format(bam_handle.filename))
        tot_mapped_reads = 0
        for line in lines:
            chrom, _len, nmapped, _nunmapped = line.split('\t')
            if chrom not in chroms_to_ignore:
                tot_mapped_reads += int(nmapped)

    else:
        tot_mapped_reads = bam_handle.mapped

    return tot_mapped_reads


def bam_blacklisted_worker(args):
    bam, chrom, start, end = args
    fh = openBam(bam)
    blacklisted = 0
    for r in fh.fetch(reference=chrom, start=start, end=end):
        if r.reference_start >= start and r.reference_start + r.infer_query_length(always=False) - 1 <= end:
            blacklisted += 1
    fh.close()
    return blacklisted


def bam_blacklisted_reads(bam_handle, chroms_to_ignore, blackListFileName=None, numberOfProcessors=1):
    blacklisted = 0
    if blackListFileName is None:
        return blacklisted

    # Get the chromosome lengths
    chromLens = {}
    lines = pysam.idxstats(bam_handle.filename)
    lines = toString(lines)
    if type(lines) is str:
        lines = lines.strip().split('\n')
    for line in lines:
        chrom, _len, nmapped, _nunmapped = line.split('\t')
        chromLens[chrom] = int(_len)

    bl = GTF(blackListFileName)
    regions = []
    for chrom in bl.chroms:
        if (not chroms_to_ignore or chrom not in chroms_to_ignore) and chrom in chromLens:
            for reg in bl.findOverlaps(chrom, 0, chromLens[chrom]):
                regions.append([bam_handle.filename, chrom, reg[0], reg[1]])

    if len(regions) > 0:
        import multiprocessing
        if len(regions) > 1 and numberOfProcessors > 1:
            pool = multiprocessing.Pool(numberOfProcessors)
            res = pool.map_async(bam_blacklisted_worker, regions).get(9999999)
        else:
            res = [bam_blacklisted_worker(x) for x in regions]
        for val in res:
            blacklisted += val

    return blacklisted
