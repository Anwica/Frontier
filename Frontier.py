#!/bin/python
import time
import numpy
import pysam
import MUSCython
import MUSCython.MultiStringBWTCython as MultiStringBWT
import os
import sys
import random
import cPickle
from multiprocessing import Pool
import argparse

# disable buffering of stdout
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

# parse commandline arguments
parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), description='Creates a frontier csv file.')
parser.add_argument('-t', '--tmp', help="path to tmp directory")
parser.add_argument('-k', '--kmerSize', type=int, default=45, help='repeated kmer size')
parser.add_argument('msbwtDir', help="path to msBWT directory")
args = parser.parse_args(sys.argv[1:])

if args.tmp is None:
    try:
        parser.tmp = os.environ["TMPDIR"]
    except KeyError: 
        print "ERROR: No tmp directory specified. Use either -t option or set TMPDIR environment variable"
        parser.print_help()
        sys.exit(1)

if args.tmp[-1] != '/':
    args.tmp += '/'

K = args.kmerSize
tmpFile = '%skmerArgs%04d.pkl' % (args.tmp, random.randint(0,9999))
    
# parallel process
def processChunk(i):
    tick = time.time()
    with open(tmpFile, 'rb') as argfile:
        filename, kmerSize, minInterval, chunk = cPickle.load(argfile)
    Huge = numpy.memmap(filename, dtype="uint8", mode="r")
    h = chunk[i]
    l = 0 if (i == 0) else chunk[i-1]
    Vec = Huge[l:h+1]
    if (Vec[0] >= kmerSize):
        delta = numpy.diff(numpy.concatenate([[0],(Vec >= kmerSize).astype('int8')]))
        lo = l + numpy.where(delta == 1)[0]
        hi = l + numpy.where(delta == -1)[0]
    else:
        delta = numpy.diff((Vec >= kmerSize).astype('int8'))
        lo = l + numpy.where(delta == 1)[0] + 1
        hi = l + numpy.where(delta == -1)[0] + 1
    intervals = numpy.vstack([lo,hi])[:,hi - lo > minInterval]
    tock = time.time()
    print "%5d of %5d, %11d:%11d %18s in %6.2f secs" % (i, len(chunk), l, h, intervals.shape, tock - tick)
    return intervals

def findKmersAtLeast(filename, kmerSize, minInterval):
    start = time.time()
    lcp = numpy.memmap(filename, dtype="uint8", mode="r")
    chunk = [i for i in sorted(random.sample(xrange(lcp.shape[0]),1000)) if lcp[i] < kmerSize]
    chunk += [lcp.shape[0]-1]
    print "Broke into %d chunks in %6.2f secs" % (len(chunk), time.time() - start)
    with open(tmpFile, 'wb') as argfile:
        cPickle.dump((filename, kmerSize, minInterval, chunk), argfile)
    start = time.time()
    job = Pool(8)
    solution = numpy.hstack(job.map(processChunk, range(len(chunk))))
    print "Done: %18s in %6.2f secs" % (solution.shape, time.time() - start)
    os.remove(tmpFile)
    return solution


lcpFilename = "%s/lcp.mmp" % args.msbwtDir

print "Scanning for repeated %d-mers" % K 
tick = time.time()
interval45 = findKmersAtLeast(lcpFilename, K, 400)
tock = time.time()
print "Finished scan. Found %s in %6.2f secs" % (interval45.shape, tock - tick)
tick = tock

msbwt = MultiStringBWT.loadBWT(args.msbwtDir)
l, h = msbwt.findIndicesOfStr('$')
print l, h, "reads"
l, h = msbwt.findIndicesOfStr('')
print l, h, "bases"
tock = time.time()
print "msBWT (%s) Loading time: %8.3f secs" % (args.msbwtDir, tock - tick)
tick = tock

lcp = numpy.memmap(lcpFilename, dtype="uint8", mode="r")
print "Opened LCP: %s %s" % (lcp.shape, lcp.dtype)

def lowComplexity(s):
    kmers = set([s[i:i+2] for i in xrange(len(s)-1)])
    return len(kmers) < 6

fp = open("%s/frontier.csv" % args.msbwtDir, "w")
fp.write("sequence,lo,hi,pre45cnt,count\n")
nondollar = msbwt.findIndicesOfStr('$')[1]
N = 0
M = 0
for index in xrange(interval45.shape[1]):
    lo = interval45[0,index]
    if (lo < nondollar):
        continue
    prefix = msbwt.recoverString(lo)[0:45]
    if (prefix.find('$') >= 0) or lowComplexity(prefix):
        continue
    hi = interval45[1,index]
    # print lo, hi, hi - lo, prefix
    Vec = lcp[lo:hi]
    delta = numpy.diff(numpy.concatenate([[0],(Vec >= 75).astype('int8'),[0]]))
    l = lo + numpy.where(delta == 1)[0]
    h = lo + numpy.where(delta == -1)[0]
    temp = numpy.logical_and(h - l > 3, h - l <= 12)
    chunk = numpy.vstack([l,h])[:,temp]
    for s, e in chunk.T:
        seq = msbwt.recoverString(s)[0:75]
        if (seq.find('$') >= 0):
            continue
        seqcount = msbwt.countOccurrencesOfSeq(MultiStringBWT.reverseComplement(seq))
        if (seqcount < 2) or (e - s + seqcount > 24):
            continue
        precount = msbwt.countOccurrencesOfSeq(MultiStringBWT.reverseComplement(prefix))
        M += 1
        fp.write("%s,%d,%d,%d,%d\n" % (seq, s, e, hi - lo + precount, e - s + seqcount))
    N += 1
    if (N % 100000 == 0):
        tock = time.time()
        print "%9d,%10d, %8.2f secs" % (N, M, tock - tick)
        tick = tock
tock = time.time()
print "Done: %9d,%10d, %8.2f secs" % (N, M, tock - tick)
fp.close()

