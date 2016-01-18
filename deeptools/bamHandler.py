import sys
import pysam


def openBam(bamFile):

    try:
        bam = pysam.Samfile(bamFile, 'rb')
    except IOError:
        sys.exit("The file {} does not exist".format(bamFile))
    except:
        sys.exit("The file {} does not have BAM format ".format(bamFile))

    try:
        if 'check_index' in dir(bam):
            assert(bam.check_index())
        else:
            # The proper check_index() function wasn't implemented until pysam 0.8.4!
            assert(bam._hasIndex())
    except:
        sys.exit("{} does not appear to have an index. You MUST index the file first!".format(bamFile))

    if bam.mapped == 0:
        sys.exit("Samtools reports that the number of mapped "
                 "reads is zero for the file {}. Please check "
                 "that the file is properly indexed and that "
                 "it contains mapped reads.".format(bamFile))

    return bam
