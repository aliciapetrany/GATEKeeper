#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from GATEKeeperUtils import GATEKeeper
import warnings

warnings.simplefilter('ignore')
parser = argparse.ArgumentParser()

#extended pipeline arguments
parser.add_argument("-f", "--fastafile", dest="file", help="Path to file containing all target sequences in fasta format", type=str)
parser.add_argument("-o", "--output", dest="output", default = "1", help="The suffix to be added to the end of all output files",type=str)
parser.add_argument("-r", "--rootpos", dest="rootpos", default=0, help="The position of the 'root' sequence, if any. Defaults to no root", type=int)
parser.add_argument("-n", "--nperms", dest="nperms", default=1000, help="The number of times to resample the minimum spanning tree", type=int)
parser.add_argument("-v", "--vcfformat", dest="vcfformat", default=0, help="Level of detail in the vcf output. 0=no vcf, 1=root vs. all variants, 2= all vs. all variants", type=int)
parser.add_argument("-m", "--metadata", dest="metadata_path", default="", help="Path to sequence metadata from the ncbi",type=str)
parser.add_argument("-x", "--testmode", dest="testmode", action="store_true", help="Run time series valdiation")
parser.add_argument("-l", "--seqlimit", dest="seqlimit", default=float("inf"), help="Cutoff on number of sequences to read in", type=int )
parser.add_argument("-t", "--threads", dest="threads", default=8, help="Number of threads", type = int)
parser.add_argument("-p", "--binpath", dest="binpath", default="bin/", help="GATEKeeper bin location", type = str)
parser.add_argument("-verb", "--verbosity", dest="verbosity", default=2, help="Verbosity switch. 0=no output, 1=warnings only, 2=warnings and info", type = int)

#Gatekeeper's arguments
parser.add_argument("-idy", "--minidentity", dest="minidentity", default=70, help="Set the minimum sequence identity (0-100) of a local alignment", type=int)
parser.add_argument("-slen", "--minslength", dest="minslength", default=15, help="Set the minimum seed length", type=int)
parser.add_argument("-alen", "--minalength", dest="minalength", default=200, help="Set the minimum alignment length", type=int)
parser.add_argument("-ind", "--maxindel", dest="maxindel", default=25, help="Set maximal indel size", type=int)
parser.add_argument("-clr", "--clustersize", dest="clustersize", default=200, help ="Set the minimum cluster size", type=int)
parser.add_argument("-u", "--unique", dest="unique", action='store_true', help="Output unique aligments only")
parser.add_argument("-sen", "--sensitive", dest="sensitive", action = "store_true", help="Sensitive mode")
parser.add_argument("-dp", "--dotplot", dest="dotplot", action="store_true", help="Output GATEKeeper dot plots")

args = parser.parse_args()

g = GATEKeeper(args.file)
g.trial_name = args.output
g.root_pos = args.rootpos
g.nperms = args.nperms
g.vcf_level = args.vcfformat
g.time_metadata_path = args.metadata_path
g.test_mode = args.testmode
g.seq_limit = args.seqlimit
g.nthreads = args.threads
g.bin_path = args.binpath
g.verbosity = args.verbosity
g.gatekeeper_args = [args.minidentity, args.minslength, args.minalength,
                     args.maxindel, args.clustersize, args.unique,
                     args.sensitive, args.dotplot]
g.run()


