import os
import sys
import pysam
import tempfile


def openBam(bamFile, bamIndex=None):

    if not os.path.exists(bamFile):
        sys.exit("Bam file {} does not exist".format(bamFile))

    if os.stat(bamFile).st_size < 10:
        sys.exit("Bam file {} is empty".format(bamFile))

    if bamIndex and bamIndex != bamFile + ".bai":
        if not os.path.exists(bamIndex):
            exit("Given Index file {} does not exists.\n"
                 "Be sure that the bam file you are using "
                 "is indexed.".format(bamIndex))

        tmpDir = tempfile.mkdtemp()
        tmpf0 = tempfile.NamedTemporaryFile(dir=tmpDir)
        tmpf0_name = tmpf0.name
        tmpf0.close()
        tmpf0bam_name = '%s.bam' % tmpf0_name
        tmpf0bambai_name = '%s.bam.bai' % tmpf0_name

        os.symlink(os.path.abspath(bamFile), tmpf0bam_name)
        os.symlink(os.path.abspath(bamIndex), tmpf0bambai_name)
        bamFile = tmpf0bam_name

    else:
        if not os.path.exists(bamFile + ".bai"):
            sys.exit("Index file {} does not exists.\n"
                     "Be sure that the bam file you are"
                     "using is indexed.".format(bamFile + ".bai"))

    try:
        bam = pysam.Samfile(bamFile, 'rb')
    except IOError:
        sys.exit("The file {} does not exits".format(bamFile))

    except:
        sys.exit("The file {} does not have BAM format ".format(bamFile))

    if bam.mapped == 0:
        sys.exit("Samtools reports that the number of mapped "
                 "reads is zero for the file {}. Please check "
                 "that the file is properly indexed and that "
                 "it contains mapped reads.".format(bamFile))

    return bam
