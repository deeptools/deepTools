import sys
import pysam


def openBam(bamFile):

    try:
        bam = pysam.Samfile(bamFile, 'rb')
    except IOError:
        sys.exit("The file {} does not exits".format(bamFile))
    except:
        sys.exit("The file {} does not have BAM format ".format(bamFile))

    try:
        assert(bam.check_index())
    except:
        sys.exit("{} does not appear to have an index. You MUST index the file first!".format(bamFile))

    if bam.mapped == 0:
        sys.exit("Samtools reports that the number of mapped "
                 "reads is zero for the file {}. Please check "
                 "that the file is properly indexed and that "
                 "it contains mapped reads.".format(bamFile))

    return bam
