import numpy as np

# own tools
from deeptools import bamHandler
from deeptools import mapReduce

old_settings = np.seterr(all='ignore')


def getFragmentLength_wrapper(args):
    return getFragmentLength_worker(*args)


def getFragmentLength_worker(chrom, start, end, bamFile, distanceBetweenBins):
    """
    Queries the reads at the given region for the distance between
    reads and the read length

    Parameters
    ----------
    chrom : str
        chromosome name
    start : int
        region start
    end : int
        region end
    bamFile : str
        BAM file name
    distanceBetweenBins : int
        the number of bases at the end of each bin to ignore

    Returns
    -------
    np.array
        an np.array, where first column is fragment length, the
        second is for read length
    """
    bam = bamHandler.openBam(bamFile)
    end = max(start + 1, end - distanceBetweenBins)
    if chrom in bam.references:
        reads = np.array([(abs(r.template_length), r.infer_query_length(always=False))
                          for r in bam.fetch(chrom, start, end)
                          if r.is_proper_pair and r.is_read1])
        if not len(reads):
            # if the previous operation produces an empty list
            # it could be that the data is not paired, then
            # we try with out filtering
            reads = np.array([(abs(r.template_length), r.infer_query_length(always=False))
                              for r in bam.fetch(chrom, start, end)])
    else:
        raise NameError("chromosome {} not found in bam file".format(chrom))

    if not len(reads):
        reads = np.array([]).reshape(0, 2)

    return reads


def get_read_and_fragment_length(bamFile, return_lengths=False, blackListFileName=None,
                                 binSize=50000, distanceBetweenBins=1000000,
                                 numberOfProcessors=None, verbose=False):
    """
    Estimates the fragment length and read length through sampling

    Parameters
    ----------
    bamFile : str
        BAM file name
    return_lengths : bool
    numberOfProcessors : int
    verbose : bool
    binSize : int
    distanceBetweenBins : int

    Returns
    -------
    d : dict
        tuple of two dictionaries, one for the fragment length and the other
        for the read length. The dictionaries summarise the mean, median etc. values
    """

    bam_handle = bamHandler.openBam(bamFile)
    chrom_sizes = list(zip(bam_handle.references, bam_handle.lengths))

    distanceBetweenBins *= 2
    fl = []
    while len(fl) < 1000 and distanceBetweenBins > 1:
        distanceBetweenBins /= 2
        stepsize = binSize + distanceBetweenBins
        imap_res = mapReduce.mapReduce((bam_handle.filename, distanceBetweenBins),
                                       getFragmentLength_wrapper,
                                       chrom_sizes,
                                       genomeChunkLength=stepsize,
                                       blackListFileName=blackListFileName,
                                       numberOfProcessors=numberOfProcessors,
                                       verbose=verbose)

        fl = np.concatenate(imap_res)

    if len(fl):
        fragment_length = fl[:, 0]
        read_length = fl[:, 1]
        if fragment_length.mean() > 0:
            fragment_len_dict = {'sample_size': len(fragment_length),
                                 'min': fragment_length.min(),
                                 'qtile25': np.percentile(fragment_length, 25),
                                 'mean': np.mean(fragment_length),
                                 'median': np.median(fragment_length),
                                 'qtile75': np.percentile(fragment_length, 75),
                                 'max': fragment_length.max(),
                                 'std': np.std(fragment_length)}
        else:
            fragment_len_dict = None

        if return_lengths and fragment_len_dict is not None:
            fragment_len_dict['lengths'] = fragment_length

        read_len_dict = {'sample_size': len(read_length),
                         'min': read_length.min(),
                         'qtile25': np.percentile(read_length, 25),
                         'mean': np.mean(read_length),
                         'median': np.median(read_length),
                         'qtile75': np.percentile(read_length, 75),
                         'max': read_length.max(),
                         'std': np.std(read_length)}
        if return_lengths:
            read_len_dict['lengths'] = read_length
    else:
        fragment_len_dict = None
        read_len_dict = None

    return fragment_len_dict, read_len_dict
