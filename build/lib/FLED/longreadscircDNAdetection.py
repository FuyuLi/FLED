import sys
import argparse
import datetime
import pysam as ps
import networkx as nx
import time
from itertools import groupby
from scipy import stats
import numpy as np
from collections import defaultdict
from progressbar import *

from FLED.ecDNAdetection import detection
from FLED.__version__ import __version__ as cm_version

class longreadscircDNAdetection :

    def __getpid__(self):

        pid = os.getpid()
        return (pid)

    def __init__(self):
        self.parser = argparse.ArgumentParser(
                description='FLED:Full Length eccDNA Detection',
            usage='''FLED <subprogram> [options]
version=%s
contact= https://github.com/FuyuLi/FLED/issues

The FLED suite

Commands:
   Detection       Identify extrachromosomal circular DNA from Nanopore reads
   Simulate        Simulate extrachromosomal circular DNA
''' % cm_version)
        subparsers = self.parser.add_subparsers()

        self.detect = subparsers.add_parser(
            name="Detection",
            description='Identify extrachromosomal circular DNA from Nanopore reads',
            prog="FLED Detection",
            usage='''FLED Detection [options]'''

        )

        if len(sys.argv) <= 1:
            self.parser.print_help()
            time.sleep(0.01)
            sys.stderr.write("\nNo argument given to FLED"
                             "\nExiting\n")
            sys.exit(0)

        else:
            if sys.argv[1] == "Detection":

                self.subprogram = self.args_detect()
                self.args = self.subprogram.parse_args(sys.argv[2:])

                object = detection(self.args.i, self.args.fq, self.args.outPrefix, self.args.directory, self.args.mapq,
                                    self.args.gap,
                                    self.args.clustering_dist, self.args.verbose,self.args.filter_level, self.args.threads,
                                    self.subprogram)
                object.run_detection()

            else:
                self.parser.print_help()
                time.sleep(0.01)
                sys.stderr.write("\nWrong argument given to FLED"
                                 "\nExiting\n")
                sys.exit(0)





    def args_detect(self):

        parser = self.detect

        parser._action_groups.pop()
        required = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')
        # prefixing the argument with -- means it's optional

        # input and output

        required.add_argument('-i', metavar='', help="Input: coordinate sorted bam file")
        required.add_argument('-fq', metavar='', help="Input: fastq file")


        if ("-i" in sys.argv) and ("-fq" in sys.argv):
            optional.add_argument('-o', '--outPrefix', metavar='',
                                  help="Prefix of Output filename",
                                  default="ecDNA_%s" % sys.argv[sys.argv.index("-i") + 1])

            optional.add_argument('-dir', '--directory', metavar='',
                                  help="Output directory, default is the working directory",
                                  default=os.getcwd())

            # mapping quality cutoff

            optional.add_argument('-q', '--mapq', type=int, metavar='',
                                           help="Minimum mapping quality allowed in the alignments. Default: 10",
                                           default=10)

            # continuous supplementary  alignments options
            optional.add_argument('-g', '--gap', type=int, metavar='',
                                  help="Maximum gap allowed in the Continuous supplementary alignments reconstruction. Default: 50",
                                  default=50)

            # soft-clipped argument
            optional.add_argument('-K', '--clustering_dist', type=int, metavar='',
                                  help="Cluster ecDNAs that are K nucleotides apart in the same interval. Default: 50",
                                  default=50)

            # verbose level

            optional.add_argument('-v', '--verbose', type=int, metavar='',
                                  help='Verbose level, 1=error,2=warning, 3=message. Default: 3',
                                  choices=[1, 2, 3], default=3)

            # filtering argument
            optional.add_argument('-F', '--filter_level', type=str, metavar='',
                                  help='Filtering level,  high=5Jreads&HighCoverage,  low=2Jreads&MiddleCoverage,  none=no-filtering. Default: low',
                                  choices=['high', 'low', 'none'], default='low')

            # run argument
            optional.add_argument('-t', '--threads', type=int, metavar='',
                                  help="Number of threads to use.Default 1", default=1)


        else:
            optional.add_argument('-o', '--outPrefix', metavar='',
                                  help="Prefix of Output filename")

            optional.add_argument('-dir', '--directory', metavar='',
                                  help="Output directory, default is the working directory",
                                  default=os.getcwd())

            # mapping quality cutoff

            optional.add_argument('-q', '--mapq', type=int, metavar='',
                                  help="Minimum mapping quality allowed in the alignments. Default: 10",
                                  default=10)

            # continuous supplementary  alignments options
            optional.add_argument('-g', '--gap', type=int, metavar='',
                                  help="Maximum gap allowed in the Continuous supplementary alignments reconstruction. Default: 50",
                                  default=50)

            # soft-clipped argument
            optional.add_argument('-K', '--clustering_dist', type=int, metavar='',
                                  help="Cluster ecDNAs that are K nucleotides apart in the same interval. Default: 50",
                                  default=50)

            # verbose level

            optional.add_argument('-v', '--verbose', type=int, metavar='',
                                  help='Verbose level, 1=error,2=warning, 3=message. Default: 3',
                                  choices=[1, 2, 3], default=3)

            # filtering argument
            optional.add_argument('-F', '--filter_level', type=str, metavar='',
                                  help='Filtering level,  high=5Jreads&HighCoverage,  low=2Jreads&MiddleCoverage,  none=no-filtering. Default: low',
                                  choices=['high', 'low', 'none'], default='low')

            # run argument
            optional.add_argument('-t', '--threads', type=int, metavar='',
                                  help="Number of threads to use.Default 1", default=1)


            parser.print_help()

            time.sleep(0.01)
            sys.stderr.write(
                "\nNo input or output input given to Detection, be sure that you are providing the flags'-i' and '-fq'"
                "\nExiting\n")
            sys.exit(0)

        # parse the commands

        if len(sys.argv[2:]) == 0:
            parser.print_help()
            time.sleep(0.01)
            sys.stderr.write("\nNo arguments given to Detection. Exiting\n")
            sys.exit(0)

        return (parser)






def main():
    run = longreadscircDNAdetection()
    pid = run.__getpid__()
    # clean
    os.system("rm -rf temp_files_%s" % pid)

if __name__ == '__main__':
    main()
