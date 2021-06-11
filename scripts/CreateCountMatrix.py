#!/bin/bash

import sys
sys.path.insert(0, "/nas/longleaf/home/anwica/packages/msbwt-0.3.0-py2.7-linux-x86_64.egg")

import pysam
import MUSCython
import MUSCython.MultiStringBWTCython as MultiStringBWT
import time
import numpy as np

def run(frontierfile, bwtfile, countfile):
	fp = open(frontierfile)
	data = fp.read().splitlines()
	fp.close()

	data = data[1:]
	k = 21
	K = 75
	m = len(data)
	n = K-k+1

	msbwt = MultiStringBWT.loadBWT(bwtfile, useMemmap=False)


	matrix = np.zeros((m,n), dtype=int)

	start = time.time()

	for i,d in enumerate(data):
		[seq,fam,c,p,s] = d.split(",")
		for j in range(0,K-k+1,1):
			probe = seq[j:j+k]
			probe_rev = MultiStringBWT.reverseComplement(probe)
			lo1, hi1 = msbwt.findIndicesOfStr(probe)
			lo2, hi2 = msbwt.findIndicesOfStr(probe_rev)
			matrix[i][j] = hi1-lo1+hi2-lo2

		if ((i+1)%50000) == 0:
			end = time.time()
			print "%d: %s" %(i,end-start)
			start = time.time()

	np.savez(countfile, count=matrix)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('frontierfile', help='List of frontier candidates')
    parser.add_argument('bwtfile', help='Path to the bwt file')
    parser.add_argument('countfile', help='Desired Path and Name of the countfile')
    
    args = parser.parse_args()
    
    frontierfile = args.frontierfile
	bwtfile = args.bwtfile
	countfile = args.countfile

    run(frontierfile, bwtfile, countfile)