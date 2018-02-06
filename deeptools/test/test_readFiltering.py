from nose.tools import assert_equal
import deeptools.estimateReadFiltering as est
import deeptools.alignmentSieve as sieve
import os.path
from os import unlink
import hashlib

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/test_data/"
BAMFILE_FILTER = ROOT + "test_filtering.bam"
BEDFILE_FILTER = ROOT + "test_filtering.blacklist.bed"
PAIREDBAMFILE_FILTER = ROOT + "test_paired.bam"


def test_estimate_read_filtering_minimal():
    """
    Minimal testing
    """
    outfile = '/tmp/test_minimal.txt'
    args = '-b {} -o {}'.format(BAMFILE_FILTER, outfile).split()
    est.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    expected = ['Sample\tTotal Reads\tMapped Reads\tAlignments in blacklisted regions\tEstimated mapped reads filtered\tBelow MAPQ\tMissing Flags\tExcluded Flags\tInternally-determined Duplicates\tMarked Duplicates\tSingletons\tWrong strand\n',
                'test_filtering.bam\t193\t193\t0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n']
    # strip the path from the output
    _ = resp[1].split("\t")
    _[0] = os.path.basename(_[0])
    resp[1] = "\t".join(_)
    assert_equal(resp, expected)
    unlink(outfile)


def test_estimate_read_filtering_params():
    """
    Minimal testing
    """
    outfile = '/tmp/test_params.txt'
    args = '-b {} --minMappingQuality 10 --samFlagExclude 512 --ignoreDuplicates -bl {} -o {}'.format(BAMFILE_FILTER, BEDFILE_FILTER, outfile).split()
    est.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()
    # strip the path from the output
    _ = resp[1].split("\t")
    _[0] = os.path.basename(_[0])
    resp[1] = "\t".join(_)
    expected = ['Sample\tTotal Reads\tMapped Reads\tAlignments in blacklisted regions\tEstimated mapped reads filtered\tBelow MAPQ\tMissing Flags\tExcluded Flags\tInternally-determined Duplicates\tMarked Duplicates\tSingletons\tWrong strand\n',
                'test_filtering.bam\t193\t193\t7\t176.1\t41.4\t0.0\t186.5\t31.6\t0.0\t0.0\t0.0\n']
    assert_equal(resp, expected)
    unlink(outfile)


def test_sieve():
    """
    Test filtering a BAM file by MAPQ, flag, and blacklist
    """
    outfile = '/tmp/test_sieve.bam'
    outfiltered = '/tmp/test_sieveFiltered.bam'
    outlog = '/tmp/test_sieve.log'
    args = '-b {} --smartLabels --minMappingQuality 10 --samFlagExclude 512 -bl {} -o {} --filterMetrics {} --filteredOutReads {}'.format(BAMFILE_FILTER, BEDFILE_FILTER, outfile, outlog, outfiltered).split()
    sieve.main(args)

    _foo = open(outlog, 'r')
    resp = _foo.readlines()
    _foo.close()

    expected = ['#bamFilterReads --filterMetrics\n',
                '#File\tReads Remaining\tTotal Initial Reads\n',
                'test_filtering\t5\t193\n']
    assert_equal(resp, expected)
    unlink(outlog)
    h = hashlib.md5(open(outfile, "rb").read()).hexdigest()
    assert(h == "977bdab227a4dbfa3fc9f27c23a3e0b7")
    unlink(outfile)

    h = hashlib.md5(open(outfiltered, "rb").read()).hexdigest()
    assert(h == "762e79b7a2245ff6b2cea4139a1455de")
    unlink(outfiltered)


def test_sieve_BED():
    """
    Test alignmentSieve with the --BED option
    """
    outfile = '/tmp/test_sieve.bed'
    args = '-b {} --minMappingQuality 10 --BED -o {}'.format(PAIREDBAMFILE_FILTER, outfile).split()
    sieve.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()

    expected = ['chr2\t5000026\t5000390\n',
                'chr2\t5000303\t5000711\n',
                'chr2\t5000384\t5000531\n',
                'chr2\t5000384\t5000531\n',
                'chr2\t5000559\t5000941\n',
                'chr2\t5000736\t5001171\n',
                'chr2\t5000819\t5001228\n',
                'chr2\t5000821\t5001158\n',
                'chr2\t5000821\t5001158\n',
                'chr2\t5000821\t5001158\n',
                'chr2\t5000834\t5001249\n',
                'chr2\t5000855\t5001277\n',
                'chr2\t5000867\t5001218\n',
                'chr2\t5000925\t5001023\n',
                'chr2\t5000925\t5001023\n',
                'chr2\t5000937\t5001338\n',
                'chr2\t5001010\t5001176\n',
                'chr2\t5001025\t5001431\n',
                'chr2\t5001050\t5001436\n',
                'chr2\t5001114\t5001413\n',
                'chr2\t5001115\t5001269\n',
                'chr2\t5001115\t5001269\n',
                'chr2\t5001226\t5001603\n',
                'chr2\t5001491\t5001527\n',
                'chr2\t5001700\t5001736\n']

    assert_equal(resp, expected)
    unlink(outfile)


def test_sieve_BED_shift():
    """
    Test alignmentSieve --BED --shift
    """
    outfile = '/tmp/test_sieve_shift.bed'
    args = '-b {} --minMappingQuality 10 --BED -o {} --shift 1 -2 3 -4'.format(PAIREDBAMFILE_FILTER, outfile).split()
    sieve.main(args)

    _foo = open(outfile, 'r')
    resp = _foo.readlines()
    _foo.close()

    expected = ['chr2\t5000027\t5000388\n',
                'chr2\t5000307\t5000708\n',
                'chr2\t5000388\t5000528\n',
                'chr2\t5000385\t5000529\n',
                'chr2\t5000560\t5000939\n',
                'chr2\t5000737\t5001169\n',
                'chr2\t5000823\t5001225\n',
                'chr2\t5000825\t5001155\n',
                'chr2\t5000825\t5001155\n',
                'chr2\t5000825\t5001155\n',
                'chr2\t5000835\t5001247\n',
                'chr2\t5000859\t5001274\n',
                'chr2\t5000868\t5001216\n',
                'chr2\t5000929\t5001020\n',
                'chr2\t5000929\t5001020\n',
                'chr2\t5000941\t5001335\n',
                'chr2\t5001011\t5001174\n',
                'chr2\t5001026\t5001429\n',
                'chr2\t5001054\t5001433\n',
                'chr2\t5001118\t5001410\n',
                'chr2\t5001119\t5001266\n',
                'chr2\t5001119\t5001266\n',
                'chr2\t5001230\t5001600\n']

    assert_equal(resp, expected)
    unlink(outfile)
