import numpy as np

# own tools
import bamHandler
import mapReduce


def getFragmentLength_wrapper(args):
    return getFragmentLength_worker(*args)


def getFragmentLength_worker(chrom, start, end, bamFile):
    bam = bamHandler.openBam(bamFile)
    end = min(end, start + 5e4)
    reads = np.array([])
    if chrom in bam.references:
        reads = np.array([abs(r.tlen)
                          for r in bam.fetch(chrom, start, end)
                          if r.is_proper_pair and r.is_read1])
    else:
        raise NameError("chromosome {} not found in bam file".format(chrom))

    return reads


def peFragmentSize(bamFile, bamFileIndex=None,
                   return_lengths=False,
                   numberOfProcessors=None, verbose=False):

    bamHandle = bamHandler.openBam(bamFile, bamFileIndex)
    chromSizes = zip(bamHandle.references, bamHandle.lengths)

    chunkSize = int(
        float(sum(bamHandle.lengths)) * 0.3 / max(numberOfProcessors,
                                                  len(bamHandle.lengths)))
    imap_res = mapReduce.mapReduce((bamHandle.filename, ),
                                   getFragmentLength_wrapper,
                                   chromSizes,
                                   genomeChunkLength=chunkSize,
                                   numberOfProcessors=numberOfProcessors,
                                   verbose=verbose)

    fl = np.concatenate(imap_res)
    if len(fl):
        fragLength = {'sample_size': len(fl),
                      'min': fl.min(),
                      'qtile25': np.percentile(fl, 25),
                      'mean': np.mean(fl),
                      'median': np.median(fl),
                      'qtile75': np.percentile(fl, 75),
                      'max': fl.max(),
                      'std': np.std(fl)}
        if return_lengths:
            fragLength['lengths'] = fl
    else:
        fragLength = None
    return fragLength
