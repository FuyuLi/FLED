import sys
import os
import argparse
import datetime
from FLED._utils import *


class detection:
    """Class for full length ecDNA detection from nanopore sequencing data"""
    def __init__(self,ont_bam,input_fq, label,out_dir,mapq_cutoff, read_gap, merge_dist,
                 verbose, filter_level, threads, parser ):
        #input-output
        self.ont_bam = ont_bam
        self.input_fq = input_fq
        self.label = label
        self.out_dir = out_dir

        #read options
        self.mapq_cutoff = mapq_cutoff
        self.read_gap = read_gap
        self.merge_dist = merge_dist

        #verbose level
        self.verbose = int(verbose)

        # filtering level
        self.filter_level = filter_level

        #parser options
        self.parser = parser

        # run options
        self.threads = threads



    def run_detection(self):

        """Function that detect full length ecDNA from ONT data """

        #OUTPUT file
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        OnesegJunction = self.out_dir + '/' + self.label + '.DiGraph.OnesegJunction.out'
        OnesegFa = self.out_dir + '/' + self.label + '.DiGraph.OnesegJunction.fa'
#        OnesegFullJunction = self.out_dir + '/' + self.label + '.DiGraph.OnesegFullJunction.out'
#        OnesegBreakJunction = self.out_dir + '/' + self.label + '.DiGraph.OnesegBreakJunction.out'
#        MulsegFullJunction = self.out_dir + '/' + self.label + '.DiGraph.MulsegFullJunction.out'
#        reads4assemblefile = self.out_dir + '/' + self.label + '.reads4canu.list'

        if self.verbose >= 3:
            print("eccDNA Detection Begins !\n")

        # timing
        begin = time.time()

        #Parse SAinfo from minimap2 aligner of circle-seq nanopore data to find continuous SA for junction detection
        SAInfo, LinearReads, readlens, readnum, mappedreadsnum = get_continuous_SAinfo(self.ont_bam,self.read_gap, self.mapq_cutoff, self.verbose, begin)

        # Full ecDNA Junction Detection
        Juncinfo, Tag, FSJ = parse_junction_from_SAinfo(SAInfo, self.merge_dist, readlens, self.verbose)
        Junctions, multiseg = JuncMerge(Juncinfo, Tag, self.merge_dist, self.verbose)

        # Full Sequence
        eccDNAseq, revisedJunc = FullSeqs(Junctions, SAInfo, self.input_fq, LinearReads, Tag, self.mapq_cutoff, self.merge_dist, self.verbose, self.threads, self.filter_level)

        # Full Single Segment Junctions Coverage
        ratiocovL, ratiocovR, mannwhitneyuL, mannwhitneyuR = Coverage_count(self.ont_bam, Junctions, revisedJunc, self.verbose, self.filter_level)  #need to multiple threads

        #OutPut
        FullJunc, BreakJunc, eccDNAnum, FulleccDNAnum = output(OnesegJunction, OnesegFa, Junctions, Tag, ratiocovL, ratiocovR, mannwhitneyuL, mannwhitneyuR, self.label, self.filter_level,
               eccDNAseq, revisedJunc, self.verbose)


        # Full Multiple Segments Junctions
        #Seginfo, Segreads = Multiple_Segment_Merge(multiseg, self.merge_dist)
        #Multiple_Segment_Write(MulsegFullJunction, Seginfo, Segreads, self.label)

        end = time.time()

        Fullread = 0
        for read in Tag :
            if (Tag[read] == 'FullBSJ') or (Tag[read] == 'FullJunc') :
                Fullread += 1


        if self.verbose >= 3:
            print("finished ecDNA Detection. Elapsed time:", (end - begin) / 60, "mins")
            print("**********Report**********")
            print(("Sample Name:        %s") % str(self.label))
            print(("Detected Reads:     %s") % str(readnum))
            print(("Mapped Reads:       %s") % str(mappedreadsnum))
            print(("Full-length Reads:  %s") % str(Fullread))
            print(("Detected eccDNA:    %s") % str(eccDNAnum))
            print(("Full-length eccDNA: %s") % str(FulleccDNAnum))
            print("**************************")
            print("Thanks for using FLED")






