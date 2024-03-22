import sys
import os
import argparse
import datetime
import subprocess
from FLED._utils import *


class detection:
    """Class for full length ecDNA detection from nanopore sequencing data"""
    def __init__(self,reffa,input_fq, label,out_dir,mapq_cutoff, read_gap, merge_dist,
                 verbose, filter_level, threads, parser ):
        #input-output
        self.reffa = reffa
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
        MulsegJunction = self.out_dir + '/' + self.label + '.DiGraph.MulsegFullJunction.out'
        MulsegFa = self.out_dir + '/' + self.label + '.DiGraph.MulsegFullJunction.fa'


        if self.verbose >= 3:
            print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"),
                    "eccDNA Detection Begins !\n")

        # timing
        begin = time.time()

        if self.verbose >= 3:
            print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"),
                    "Minimap2 alignment Begins !\n")

        if not os.path.exists(self.out_dir + '/MappingResult'):
            os.makedirs(self.out_dir + '/MappingResult')
        samfile = self.out_dir + '/MappingResult/' + self.label + ".sam"
        bamfile = self.out_dir + '/MappingResult/' + self.label + ".sorted.bam"
        baifile = self.out_dir + '/MappingResult/' + self.label + ".sorted.bam.bai"

        alignment = subprocess.call(["minimap2", "-t", str(self.threads), "-ax", "map-ont", self.reffa, self.input_fq, "-o", samfile], shell=False)
        if alignment  != 0:
            print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"),
                  "An error happened during minimap2 for reads alignment. Exiting")
            sys.exit()
        else:
            print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"),
                    "Minimap2 alignment Done!")
        samtoolsSort = subprocess.call(["samtools", "sort", "-@", str(self.threads), "-O", "bam", "-o", bamfile, samfile], shell=False)
        samtoolsIndex = subprocess.call(["samtools", "index",bamfile, baifile], shell=False)
        if (samtoolsSort + samtoolsIndex) != 0:
            print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"),
                  "An error happened during samtools. Exiting")
            sys.exit()


        self.ont_bam = bamfile

        #Parse SAinfo from minimap2 aligner of circle-seq nanopore data to find continuous SA for junction detection
        SAInfo, LinearReads, readlens, readnum, mappedreadsnum = get_continuous_SAinfo(self.ont_bam,self.read_gap, self.mapq_cutoff, self.verbose, begin)

        # Full ecDNA Junction Detection
        Juncinfo, Tag, FSJ = parse_junction_from_SAinfo(SAInfo, self.merge_dist, readlens, self.verbose)
        Junctions, multiseg = JuncMerge(Juncinfo, Tag, self.merge_dist, self.verbose)

        
        if len(Junctions) > 0 :
            # Full Sequence
            eccDNAseq, revisedJunc = FullSeqs(Junctions, SAInfo, self.input_fq, LinearReads, Tag, self.mapq_cutoff, self.merge_dist, self.verbose, self.threads, self.filter_level)

            # Full Single Segment Junctions Coverage
            covfile = self.out_dir + '/MappingResult/' + self.label + ".coverage"

            samtoolsIndex = subprocess.call(["samtools", "depth", "-Q", str(self.mapq_cutoff), "-d", "0", "-o", covfile, bamfile], shell=False)
            if samtoolsIndex  != 0:
                print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"),
                  "An error happened during samtools depth. Exiting")
                sys.exit()

            ratiocovL, ratiocovR, mannwhitneyuL, mannwhitneyuR = Coverage_count(covfile, Junctions, revisedJunc, self.verbose, self.filter_level)  #need to multiple threads

            #OutPut
            FullJunc, BreakJunc, eccDNAnum, FulleccDNAnum = output(OnesegJunction, OnesegFa, Junctions, Tag, ratiocovL, ratiocovR, mannwhitneyuL, mannwhitneyuR, self.label, self.filter_level, eccDNAseq, revisedJunc, self.verbose)
        else :
            eccDNAnum = 0
            FulleccDNAnum = 0


        if len(multiseg) > 0 :
            # Full Multiple Segments Junctions
            Seginfo, Segreads = Multiple_Segment_Merge(multiseg, self.merge_dist)

            # Map full-length reads to PseudoReference for Multiple Segments Junctions
            tmpbam, tmpfq = MS_PseudoReference(self.out_dir, Segreads, self.reffa, self.input_fq, self.threads, self.verbose)


            multibegin = time.time()
            MS_SAInfo, MSLinearReads, MSreadlens, MSreadnum, MSmappedreadsnum = get_continuous_SAinfo(tmpbam, self.read_gap, self.mapq_cutoff, self.verbose, multibegin)

            # Multiple-segments eccDNA detecion
            MS_Juncinfo, MS_Junctions, MS_Tag = parse_junction_from_SAinfo_of_MultipleSegment_FLreads(MS_SAInfo, MSreadlens, self.verbose)

            # Full Sequence
            MSeccDNAseq = FullSeqs_MS(MS_Junctions, MS_SAInfo, tmpfq, self.verbose, self.filter_level)

            # OutPut
            FullMSeccDNAnum = MSoutput(MulsegJunction, MulsegFa, MS_Junctions, MSeccDNAseq, self.label, self.filter_level, self.verbose)
        else :
            FullMSeccDNAnum = 0

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
            print(("Full-length multiple-segments eccDNA: %s") % str(FullMSeccDNAnum))
            print("**************************")
            print("Thanks for using FLED")






