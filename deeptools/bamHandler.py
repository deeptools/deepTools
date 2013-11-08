import pysam
import tempfile, os

def openBam(bamFile, bamIndex=None):
    try:
        open(bamFile, 'r').close()
    except IOError:
        print "Bam file {} does not exist".format(bamFile)
        exit()
    if bamIndex and bamIndex != bamFile + ".bai":
        try:
            f = open(bamIndex, 'r')
            bamIndex = f.name
            f.close()
        except IOError:
            print "Given Index file {} does not exists".format(bamIndex)
            print "Be sure that the bam file you are using is indexed."
            exit()

        tmpDir = tempfile.mkdtemp()
        tmpf0 = tempfile.NamedTemporaryFile( dir=tmpDir )
        tmpf0_name = tmpf0.name
        tmpf0.close()
        tmpf0bam_name = '%s.bam' % tmpf0_name
        tmpf0bambai_name = '%s.bam.bai' % tmpf0_name
        
        os.symlink( bamFile, tmpf0bam_name )
        os.symlink( bamIndex, tmpf0bambai_name )
        bamFile = tmpf0bam_name

    else:
        try:
            open(bamFile+".bai", 'r').close()
        except IOError:
            print "Index file {} does not exists".format(bamFile+".bai")
            print "Be sure that the bam file you are using is indexed."
            exit()
        
    try:
        bam = pysam.Samfile(bamFile, 'rb')
    except IOError:
        print "The file {} does not exits".format(bamFile)
        exit()

    except:
        print "The file {} does not have BAM format ".format(bamFile)
        exit()
        
         
    return bam
