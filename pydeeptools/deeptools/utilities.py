import sys
import os
from deeptoolsintervals import GTF
from deeptools.bamHandler import openBam
import matplotlib as mpl
mpl.use('Agg')
import numpy as np


debug = 0


def smartLabel(label):
    """
    Given a file name, likely with a path, return the file name without the path
    and with the file extension removed. Thus, something like /path/to/some.special.file
    should return some.special, since only the first extension (if present)
    should be stripped.
    """
    lab = os.path.splitext(os.path.basename(label))[0]
    if lab == '':
        # Maybe we have a dot file?
        lab = os.path.basename(label)
    return lab


def smartLabels(labels):
    return [smartLabel(x) for x in labels]


def convertCmap(c, vmin=0, vmax=1):
    cmap = mpl.cm.get_cmap(c)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap_rgb = []

    for i in range(255):
        k = mpl.colors.colorConverter.to_rgb(cmap(norm(i)))
        cmap_rgb.append(k)

    h = 1.0 / 254
    colorScale = []
    for k in range(255):
        C = list(map(np.uint8, np.array(cmap(k * h)[:3]) * 255))
        colorScale.append([k * h, 'rgb' + str((C[0], C[1], C[2]))])

    return colorScale


def getTLen(read, notAbs=False):
    """
    Get the observed template length of a read. For a paired-end read, this is
    normally just the TLEN field. For SE reads this is the observed coverage of
    the genome (excluding splicing).
    """
    if abs(read.template_length) > 0:
        if notAbs:
            return read.template_length
        return abs(read.template_length)

    tlen = 0
    try:
        # the cigartuples property apparently didn't always exist
        for op, opLen in read.cigartuples:
            if op == 0:
                tlen += opLen
            elif op == 2:
                tlen += opLen
            elif op == 7:
                tlen += opLen
            elif op == 8:
                tlen += opLen
    except:
        pass

    return tlen


def getCommonChrNames(bamFileHandles, verbose=True):
    r"""
    Compares the names and lengths of a list of bam file handles.
    The input is list of pysam file handles.

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
        try:
            # BAM file
            return [(x, y) for x, y in zip(bam_handler.references, bam_handler.lengths)]
        except:
            return [(k, v) for k, v in bam_handler.chroms().items()]

    def print_chr_names_and_size(chr_set):
        sys.stderr.write("chromosome\tlength\n")
        for name, size in chr_set:
            sys.stderr.write("{0:>15}\t{1:>10}\n".format(name, size))

    common_chr = set(get_chrom_and_size(bamFileHandles[0]))
    non_common_chr = set()

    for j in range(1, len(bamFileHandles)):
        _names_and_size = set(get_chrom_and_size(bamFileHandles[j]))
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
                                 "lengths from file\n{}\n".format(bamFileHandles.name))
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
    for tuple in get_chrom_and_size(bamFileHandles[0]):
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
    Return a temporary file name. The calling function is responsible for
    deleting this upon completion.
    """
    import tempfile
    _tempFile = tempfile.NamedTemporaryFile(prefix="_deeptools_",
                                            suffix=suffix,
                                            delete=False)

    memFileName = _tempFile.name
    _tempFile.close()
    return memFileName


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


def bam_total_reads(bam_handle, chroms_to_ignore, stats):
    """
    Count the total number of mapped reads in a BAM file, filtering
    the chromosome given in chroms_to_ignore list
    """
    if chroms_to_ignore:
        return sum([s[0] for k, s in stats.items() if k not in chroms_to_ignore])
    else:
        return sum([s[0] for s in stats.values()])


def bam_blacklisted_worker(args):
    bam, chrom, start, end = args
    fh = openBam(bam)
    blacklisted = 0
    for r in fh.fetch(reference=chrom, start=start, end=end):
        if r.is_unmapped:
            continue
        if r.reference_start >= start and r.reference_start + r.infer_query_length(always=False) - 1 <= end:
            blacklisted += 1
    fh.close()
    return blacklisted


def bam_blacklisted_reads(bam_handle, chroms_to_ignore, blackListFileName=None, numberOfProcessors=1):
    blacklisted = 0
    if blackListFileName is None:
        return blacklisted

    # Get the chromosome lengths
    chromLens = {x: y for x, y in zip(bam_handle.references, bam_handle.lengths)}

    bl = GTF(blackListFileName)
    hasOverlaps, minOverlap = bl.hasOverlaps(returnDistance=True)
    if hasOverlaps:
        sys.exit("Your blacklist file(s) has (have) regions that overlap. Proceeding with such a file would result in deepTools incorrectly calculating scaling factors. As such, you MUST fix this issue before being able to proceed.\n")
    if minOverlap < 1000:
        sys.stderr.write("WARNING: The minimum distance between intervals in your blacklist is {}. It makes little biological sense to include small regions between two blacklisted regions. Instead, these should likely be blacklisted as well.\n".format(minOverlap))

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
