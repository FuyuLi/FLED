import pysam as ps
import networkx as nx
import time
from itertools import groupby
import datetime
from scipy import stats
import numpy as np
from collections import defaultdict
from progressbar import *
from Bio import SeqIO
from spoa import poa
import sys
import argparse
from tqdm import tqdm
import subprocess
from Bio.SeqRecord import SeqRecord
from queue import Queue


class Segment(object):
    '''
    Modified from https://github.com/brentp/cigar
    '''
    def __init__(self, pos, cigar):
        self.ref_start = int(pos) - 1
        self.ref_end = self.ref_start
        read_consuming_ops = ("M", "I")
        ref_consuming_ops = ("M", "D")
        cig_iter = groupby(cigar, lambda c: c.isdigit())
        self.read_start, self.read_end = 0, 0
        for i, (g, n) in enumerate(cig_iter):
            counts, tag = int("".join(n)), "".join(next(cig_iter)[1])
            if i == 0 and tag == 'S':
                self.read_start += counts
                self.read_end += counts
            if tag in read_consuming_ops:
                self.read_end += counts
            if tag in ref_consuming_ops:
                self.ref_end += counts

def longest_path(G):
    dist = {} # stores [node, distance] pair
    for node in nx.topological_sort(G):
        # pairs of dist,node for all incoming edges
        pairs = [(dist[v][0]+1,v) for v in G.pred[node]]
        if pairs:
            dist[node] = max(pairs)
        else:
            dist[node] = (0, node)
    node,(length,_)  = max(dist.items(), key=lambda x:x[1])
    path = []
    while length > 0:
        path.append(node)
        length,node = dist[node]
    return list(reversed(path))

def get_continuous_SAinfo(ont_bam,read_gap, mapq_cutoff, verbose, begin):
    '''
    Parse SAinfo from minimap2 aligner of circle-seq nanopore data to find continuous SA for junction detection
    '''
    if verbose >= 3:
        print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Parse SAinfo for continuous SA reads\n")
    bamFile = ps.AlignmentFile("%s" % ont_bam, "rb")
    valid_cutoff = 0.7
    SAInfo = {}
    LinearReads = defaultdict(dict)
    readlens = {}
    processed_reads = 0
    readlist = {}
    for read in bamFile:
        readid = read.query_name
        if readid not in readlist.keys() :
            processed_reads += 1
            readlist[readid] = processed_reads
        if verbose >= 3:
            if (processed_reads / 100000).is_integer() == True:
                partial_timer = time.time()
                partial_time = (partial_timer - begin) / 60
                print("Processed %s reads in %s mins" % (processed_reads, round(partial_time, 3)))
                processed_reads += 1
        if read.is_unmapped:  # unmapped reads
            continue
        if read.is_supplementary:  # supplementary reads
            continue
        if read.mapq < mapq_cutoff:
            continue
        chr1 = bamFile.get_reference_name(read.reference_id)
        strand1 = '+' if not read.is_reverse else '-'
        readlen = read.infer_read_length()
        readlens[readid] = readlen
        if not read.has_tag('SA'):
            loc = [read.query_alignment_start, read.query_alignment_end, chr1, strand1,
                   read.reference_start, read.reference_end, read.query_alignment_length]
            alignedlen = read.query_alignment_length
            if (alignedlen/readlen >= valid_cutoff):
                LinearReads[chr1][readid] = loc
        else:
            saInfo = read.get_tag('SA').split(';')[:-1]
            loc = [read.query_alignment_start, read.query_alignment_end, chr1, strand1,
                   read.reference_start, read.reference_end]
            segments = [loc]
            for sa in saInfo:
                chr2, pos, strand2, cigar = sa.split(',')[:4]
                segment = Segment(pos=pos, cigar=cigar)
                segments.append([segment.read_start, segment.read_end, chr2, strand2,
                                 segment.ref_start, segment.ref_end])
            segments.sort()
            SAgraph = nx.DiGraph()
            SAgraph.add_node('END')
            nodelist = {}
            for item in segments:
                nodeID = '-'.join([str(item[0]), str(item[1])])
                nodelist[nodeID] = item
                SAgraph.add_node(nodeID, segS=item[0], segE=item[1], chrom=item[2], strand=item[3], refS=item[4],
                                 refE=item[5], length=item[1] - item[0])
                SAgraph.add_edge(nodeID, 'END', weight=(item[1] - item[0]))
            checkedseg = []
            for item in segments:
                for item2 in checkedseg:
                    if (abs(item2[1] - item[0]) <= read_gap) and item[1] > item2[1]:
                        node1 = '-'.join([str(item[0]), str(item[1])])
                        node2 = '-'.join([str(item2[0]), str(item2[1])])
                        SAgraph.add_edge(node2, node1, weight=(item2[1] - item2[0]))
                checkedseg.append(item)
            LPath = nx.dag_longest_path(SAgraph, weight='weight')
            if len(LPath) == 2:
                linearloc = nodelist[LPath[0]]
                alignedlen = linearloc[1] - linearloc[0]
                if (alignedlen / readlen >= valid_cutoff):
                    linearloc.append(alignedlen)
                    LinearReads[chr1][readid] = linearloc
            elif len(LPath) > 2:
                info = []
                alignedlen = 0
                for item in LPath:
                    if item != 'END':
                        alignedlen += SAgraph.nodes[item]['length']
                        info.append(
                            [SAgraph.nodes[item]['segS'], SAgraph.nodes[item]['segE'], SAgraph.nodes[item]['chrom'],
                             SAgraph.nodes[item]['strand'],
                             SAgraph.nodes[item]['refS'], SAgraph.nodes[item]['refE'], SAgraph.nodes[item]['length']])
                if ((alignedlen / readlen) >= valid_cutoff) :
                    SAInfo[readid] = info
    readnum = len(readlist)
    mappedreadsnum = readnum - bamFile.unmapped
    return (SAInfo, LinearReads, readlens, readnum, mappedreadsnum)

#SAInfo, LinearReads, readlens, readnum, mappedreadsnum = get_continuous_SAinfo(ont_bam,read_gap, mapq_cutoff, verbose, begin)


class AlignSeg(object):
    '''
    aligned segment information : [readS, readE, chrom, strand, refS, refE, length]
    '''
    def __init__(self, info):
        self.readS = info[0]
        self.readE = info[1]
        self.chrom = info[2]
        self.strand = info[3]
        self.refS = info[4]
        self.refE = info[5]
        self.len = info[6]

def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        >>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        >>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    data.sort()
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups

def parse_junction_from_SAinfo(SAInfo, merge_dist, readlens, verbose) :
    '''
    Parse junction from continuous SA for junction detection to Tag_Full & Tag_ToAssembly(BreakBSJ + FSJ + interchrom + Unclassified)
    '''
    if verbose >= 3:
        print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Parse continuous SAinfo for junction detection\n")
    Tag = {}    # Full, Break, FullJunc, FSJ, invalid(Interchrom, Diffstrand, Unclassified)
    # FullBSJ, BreakBSJ, FSJ, Interchrom, Unclassified, Diffstrand, Mixed
    '''
    FullBSJ        # Reads with Full sequence between BSJ
    BreakBSJ       # 2 SA Reads with BSJ but without full sequence
    FSJ            # 2 SA Reads with FSJ
    Interchrom     # 2 SA reads with interchrom junction -> incorrect junction by alignment
    Unclassified   # 2 SA reads with unknown situation
    Diffstrand     # 2 SA reads with same chrom but different strand -> incorrect junction by alignment
    Mixed          # >=3 SA reads without BSJ
    FullJunc       # >=3 SA reads with Full junction
    '''
    Juncinfo = {}
    FSJ = {}
    for readid in SAInfo.keys() :
        if len(SAInfo[readid]) < 2 :
            print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"),"An error happenend during execution. Exiting")
            sys.exit()
        elif len(SAInfo[readid]) == 2 :
            segment1 = AlignSeg(SAInfo[readid][0])
            segment2 = AlignSeg(SAInfo[readid][1])
            if segment1.chrom != segment2.chrom :
                Tag[readid] = 'invalid'
            else :
                if segment1.strand != segment2.strand :
                    Tag[readid] = 'invalid'
                else :
                    if segment1.refS > segment2.refE :
                        Tag[readid] = 'BreakBSJ'
                        Juncinfo[readid] = [','.join(['S', segment2.chrom, segment2.strand, str(segment2.refS)]), ','.join(['E', segment1.chrom, segment1.strand, str(segment1.refE)])]
                    elif (segment1.refE >= segment2.refE) and (segment2.refS <= segment1.refS <= segment2.refE) :
                        Juncinfo[readid] = [','.join(['S', segment2.chrom, segment2.strand, str(segment2.refS)]), ','.join(['E', segment1.chrom, segment1.strand, str(segment1.refE)])]
                        Tag[readid] = 'FullBSJ'
                    elif segment1.refE < segment2.refS :
                        FSJ[readid] = ','.join([segment2.chrom, (str(segment1.refS) + '-' + str(segment1.refE)),  (str(segment2.refS) + '-' + str(segment2.refE)), segment2.strand])
                        Tag[readid] = 'FSJ'
                    else :
                        Tag[readid] = 'invalid'
        else :
            Start = {}
            End = {}
            nodeS = {}
            nodeE = {}
            info = SAInfo[readid]
            SAnum = len(info)
            End[','.join([info[0][2],info[0][3]])] = [info[0][5]]
            Start[','.join([info[SAnum-1][2],info[SAnum-1][3]])] = [info[SAnum-1][4]]
            for subscript in range(1, SAnum - 1):
                seginfo = info[subscript]
                if ','.join([seginfo[2],seginfo[3]]) in Start :
                    Start[','.join([seginfo[2],seginfo[3]])].append(seginfo[4])
                else :
                    Start[','.join([seginfo[2], seginfo[3]])] = [seginfo[4]]
                if ','.join([seginfo[2],seginfo[3]]) in End :
                    End[','.join([seginfo[2],seginfo[3]])].append(seginfo[5])
                else :
                    End[','.join([seginfo[2], seginfo[3]])] = [seginfo[5]]
            for chrom in Start :
                clusters = cluster(Start[chrom],merge_dist)
                for clusteri in clusters :
                    modei = stats.mode(clusteri)[0][0]
                    for starti in clusteri :
                        nodeS[','.join([chrom,str(starti)])] = ','.join(['S',chrom,str(modei)])
            for chrom in End :
                clusters = cluster(End[chrom],merge_dist)
                for clusteri in clusters :
                    modei = stats.mode(clusteri)[0][0]
                    for starti in clusteri :
                        nodeE[','.join([chrom,str(starti)])] = ','.join(['E',chrom,str(modei)])
            DG = nx.DiGraph()
            lastnode = nodeE[','.join([info[0][2],info[0][3],str(info[0][5])])]
            for subscript in range(1, SAnum - 1):
                seginfo = info[subscript]
                nodestart = nodeS[','.join([seginfo[2],seginfo[3],str(seginfo[4])])]
                nodeend = nodeE[','.join([seginfo[2],seginfo[3],str(seginfo[5])])]
                if DG.has_edge(lastnode, nodestart) :
                    DG.edges[lastnode, nodestart]['weight'] += 1
                else :
                    DG.add_edge(lastnode,nodestart,weight=1)
                DG.add_edge(nodestart,nodeend,weight=1000)
                lastnode = nodeend
            seginfo = info[SAnum - 1]
            nodestart = nodeS[','.join([seginfo[2], seginfo[3], str(seginfo[4])])]
            if DG.has_edge(lastnode, nodestart):
                DG.edges[lastnode, nodestart]['weight'] += 1
            else :
                DG.add_edge(lastnode, nodestart, weight=1)
            cycles = list(nx.simple_cycles(DG))
            maxweight = 0
            for cycle in cycles :
                if len(cycle)%2 == 0 :
                    segnum = len(cycle)/2
                    cycle.append(cycle[0])
                    distance = []
                    for i in range(len(cycle) - 1):
                        if 'E,' in cycle[i] :
                            get_nodes_weight = DG.get_edge_data(cycle[i], cycle[i + 1])
                            weight = get_nodes_weight['weight']
                            distance.append(weight)
                    cycleweight = sum(distance)/segnum
                    if cycleweight > maxweight :
                        maxweight = cycleweight
                        bestcycle = cycle
            if maxweight == 0 :
                Tag[readid] = 'invalid'
            else :
                cyclelen = 0
                for i in range(len(bestcycle) - 1):
                    if 'S,' in bestcycle[i]:
                        if 'E,' in bestcycle[i + 1] :
                            if (bestcycle[i].split(',')[1] == bestcycle[i+1].split(',')[1]) and (bestcycle[i].split(',')[2] == bestcycle[i+1].split(',')[2]) :
                                cyclelen += (int(bestcycle[i+1].split(',')[3])-int(bestcycle[i].split(',')[3]))
                            else :
                                Tag[readid] = 'FullJunc'
                                cyclelen = 0
                                break
                cyclelen = cyclelen * maxweight
                if cyclelen/readlens[readid] >= 0.5 :
                    bestcycle.pop()
                    Juncinfo[readid] = bestcycle
                    Tag[readid] = 'FullJunc'
                else :
                    Tag[readid] = 'invalid'
    return (Juncinfo, Tag, FSJ)

#Juncinfo, Tag, FSJ = parse_junction_from_SAinfo(SAInfo, merge_dist, readlens, verbose)


def segments_cluster(data, maxgap):
    '''
    Arrange data into groups where successive elements
    differ by no more than *maxgap*
    '''
    data.sort()
    groups = [[data[0]]]
    lastchr = data[0][0]
    chrnum = 0
    for x in range(1,len(data)) :
        segment = data[x]
        if not segment[0] == lastchr :
            chrnum = len(groups)
            lastchr = segment[0]
        newGroup = True
        for group in groups[chrnum:] :
            if (segment[0] == group[-1][0]) and (abs(segment[1] - group[-1][1]) <= maxgap) and (abs(segment[2] - group[-1][2]) <= maxgap):
                group.append(segment)
                newGroup = False
                break
        if newGroup:
            groups.append([segment])
    return groups



def JuncMerge(Juncinfo, Tag, merge_dist, verbose) :
    '''
    :param Juncinfo:
    :param SAInfo:
    :param Tag:
    :return:
    No need to take strand information into account while merging
    '''
    if verbose >= 3 :
        print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Merge and Count junctions\n")
    onesegments = []
    multiseg = defaultdict(list)
    for readid in Juncinfo :
        if len(Juncinfo[readid]) == 2 :
            if Juncinfo[readid][0].split(',')[0] == 'S' :
                segment = [Juncinfo[readid][0].split(',')[1],int(Juncinfo[readid][0].split(',')[3]),
                           int(Juncinfo[readid][1].split(',')[3]), readid, Juncinfo[readid][0].split(',')[2]]
            elif Juncinfo[readid][0].split(',')[0] == 'E' :
                segment = [Juncinfo[readid][1].split(',')[1],int(Juncinfo[readid][1].split(',')[3]),
                           int(Juncinfo[readid][0].split(',')[3]), readid, Juncinfo[readid][1].split(',')[2]]
            else :
                Tag[readid] = 'invalid'
                print(readid)
                continue
            onesegments.append(segment)
        elif len(Juncinfo[readid])%2 == 0 :
            segmentNum = int(len(Juncinfo[readid])/2)
            segments = []
            if Juncinfo[readid][0].split(',')[0] == 'S':
                for st in range(0,segmentNum) :
                    segmenti = []
                    segmenti.append(Juncinfo[readid][2*st].split(',')[1])
                    segmenti.append(Juncinfo[readid][2 * st].split(',')[3])
                    segmenti.append(Juncinfo[readid][2 * st + 1].split(',')[3])
                    segmenti.append(Juncinfo[readid][2 * st + 1].split(',')[2])
                    segments.append(segmenti)
            elif Juncinfo[readid][0].split(',')[0] == 'E':
                if (Juncinfo[readid][0].split(',')[1] == Juncinfo[readid][-1].split(',')[1]) and (Juncinfo[readid][0].split(',')[2] == Juncinfo[readid][-1].split(',')[2]) :
                    segmenti = []
                    segmenti.append(Juncinfo[readid][-1].split(',')[1])
                    segmenti.append(Juncinfo[readid][-1].split(',')[3])
                    segmenti.append(Juncinfo[readid][0].split(',')[3])
                    segmenti.append(Juncinfo[readid][0].split(',')[2])
                    segments.append(segmenti)
                    for st in range(0, segmentNum-1):
                        segmenti = []
                        segmenti.append(Juncinfo[readid][2 * st + 1].split(',')[1])
                        segmenti.append(Juncinfo[readid][2 * st + 1].split(',')[3])
                        segmenti.append(Juncinfo[readid][2 * st + 2].split(',')[3])
                        segmenti.append(Juncinfo[readid][2 * st + 2].split(',')[2])
                        segments.append(segmenti)
            else:
                Tag[readid] = 'invalid'
                continue
            if segmentNum == 2 :
                segments.sort()
            segments.append([readid])
            multiseg[segmentNum].append(segments)
        else :
            Tag[readid] = 'invalid'
            continue
    onesegGroup_iteration1 = segments_cluster(onesegments, merge_dist)
    Junctions = {}
    for clusteri in onesegGroup_iteration1 :
        chrom = clusteri[0][0]
        starts = []
        ends = []
        reads = []
        strands = []
        for seg in clusteri :
            starts.append(seg[1])
            ends.append(seg[2])
            reads.append(seg[3])
            strands.append(seg[4])
        StartPos = stats.mode(starts)[0][0]
        EndPos = stats.mode(ends)[0][0]
        Strand = stats.mode(strands)[0][0]
        junc = '\t'.join([chrom,  str(StartPos), str(EndPos), Strand])
        Junctions[junc] = reads
    return(Junctions, multiseg)

#Junctions, multiseg = JuncMerge(Juncinfo, Tag, merge_dist, verbose)



def rev(seq):
    base_trans = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N', 'a':'t', 'c':'g', 't':'a', 'g':'c', 'n':'n'}
    rev_seq = list(reversed(seq))
    rev_seq_list = [base_trans[k] for k in rev_seq]
    rev_seq = ''.join(rev_seq_list)
    return(rev_seq)

class JuncInfo(object):
    '''
    aligned segment information : [chrom, startpos, endpos, strand]
    '''
    def __init__(self, info):
        self.chrom = info[0]
        self.startpos = int(info[1])
        self.endpos = int(info[2])
        self.strand = info[3]
        self.len = self.endpos - self.startpos

def subseq(Spos, Epos, Strand, cutoff, junc, segseqs, readseq, readid):
    if (Epos - Spos) / junc.len >= cutoff:
        if Strand == '+' :
            seq = readseq[readid][Spos:Epos]
        else :
            seq = rev(readseq[readid])[Spos:Epos]
        segseqs.append(seq)
    return(segseqs)


def split_bins(tasks, nth):
    '''
    split the size of data into sections for multi-processes
    :param tasks: [list] total tasks of multi-processes
    :param nth: [int] number of threads
    :return: sub_tasks [list in list] the amount of tasks performed by each child process
     
    '''
    length, tasks = len(tasks), list(tasks)
    bin_size, sub_tasks = int(length / nth), []
    if length % nth: bin_size += 1
    for idx in range(nth):
        sub_tasks.append([])
    for task in range(length) :
        subscript = task%nth
        sub_tasks[subscript].append(tasks[task])
    return sub_tasks


def consensusSeq(tasks, kargs):
    if kargs['filtering_level'] == 'high' :
        supporting_reads = 5
    elif kargs['filtering_level'] == 'low' :
        supporting_reads = 2
    else :
        supporting_reads = 1
    revisedJunc = {}
    eccDNAseq = defaultdict(str)
    for junction in tqdm(tasks, desc=' pid:' + str(os.getpid())):
        FullJuncNum = 0
        FullBSJNum = 0
        BreakBSJNum = 0
        for read in kargs['Junctions'][junction]:
            if kargs['Tag'][read] == 'FullJunc':
                FullJuncNum += 1
            elif kargs['Tag'][read] == 'FullBSJ':
                FullBSJNum += 1
            elif kargs['Tag'][read] == 'BreakBSJ':
                BreakBSJNum += 1
        if FullBSJNum == 0 and FullJuncNum == 0:  #only BreakBSJ
            ccs = '-'
            revisedJunc[junction] = junction
            eccDNAseq[junction] = ccs
        elif (len(kargs['Junctions'][junction]) < supporting_reads):
            ccs = '-'
            revisedJunc[junction] = junction
            eccDNAseq[junction] = ccs
        else :
            segseqs = []
            revisedend = []
            junc = JuncInfo(junction.split('\t'))
            for readid in kargs['LinearReads'][junc.chrom].keys():
                segment = AlignSeg(kargs['LinearReads'][junc.chrom][readid])
                if (junc.chrom == segment.chrom) and (junc.startpos <= segment.refS) and (junc.endpos >= segment.refE):
                    segseqs = subseq(segment.readS, segment.readE, segment.strand, 0.75, junc, segseqs, kargs['readseq'], readid)
            for readid in kargs['Junctions'][junction] :
                sainfo = kargs['SAInfo'][readid]
                if kargs['Tag'][readid] == 'FullBSJ' :
                    segment = AlignSeg(sainfo[0])
                    nextseg = AlignSeg(sainfo[1])
                    misbase = nextseg.readS - segment.readE
                    if misbase >= 0:
                        revisedend.append(junc.endpos)
                    else:
                        revisedend.append((junc.endpos + misbase))
                    segseqs = subseq(segment.readS, nextseg.readS, segment.strand, 0.2, junc, segseqs, kargs['readseq'], readid)
                    segseqs = subseq(nextseg.readS, nextseg.readE, segment.strand, 0.2, junc, segseqs, kargs['readseq'], readid)
                elif kargs['Tag'][readid] == 'BreakBSJ' :
                    segment = AlignSeg(sainfo[0])
                    nextseg = AlignSeg(sainfo[1])
                    misbase = nextseg.readS - segment.readE
                    if misbase >= 0:
                        revisedend.append(junc.endpos)
                    else:
                        revisedend.append((junc.endpos + misbase))
                    segseqs = subseq(segment.readS, nextseg.readS, segment.strand, 0.3, junc, segseqs, kargs['readseq'], readid)
                    segseqs = subseq(nextseg.readS, nextseg.readE, segment.strand, 0.3, junc, segseqs, kargs['readseq'], readid)
                elif kargs['Tag'][readid] == 'FullJunc' :
                    for sanum in range(len(sainfo) - 1):
                        segment = AlignSeg(sainfo[sanum])
                        nextseg = AlignSeg(sainfo[sanum + 1])
                        if (junc.chrom == segment.chrom) and (junc.startpos - kargs['merge_dist'] <= segment.refS) and (
                                    abs(junc.endpos - segment.refE) <= kargs['merge_dist']):
                            if (junc.chrom == nextseg.chrom) and (abs(junc.startpos - nextseg.refS) <= kargs['merge_dist']) and (
                                        junc.endpos + kargs['merge_dist'] >= nextseg.refE):
                                misbase = nextseg.readS - segment.readE
                                if misbase >= 0:
                                    revisedend.append(junc.endpos)
                                else:
                                    revisedend.append((junc.endpos) + misbase)
                                segseqs = subseq(segment.readS, nextseg.readS, segment.strand, 0.5, junc, segseqs, kargs['readseq'], readid)
                            elif (abs(junc.startpos - segment.refS) <= kargs['merge_dist']) :
                                segseqs = subseq(segment.readS, segment.readE, segment.strand, 0.75, junc, segseqs, kargs['readseq'], readid)
                    sanum += 1
                    segment = AlignSeg(sainfo[sanum])
                    if (junc.chrom == segment.chrom) and (junc.endpos + kargs['merge_dist'] >= segment.refE) and (
                                abs(junc.startpos - segment.refS) <= kargs['merge_dist']):
                        segseqs = subseq(segment.readS, segment.readE, segment.strand, 0.5, junc, segseqs, kargs['readseq'], readid)
            if len(revisedend) > 0 :
                newend = stats.mode(revisedend)[0][0]
            else :
                newend = junc.endpos
            newjunc = '\t'.join([junc.chrom, str(junc.startpos), str(newend), junc.strand])
            revisedJunc[junction] = newjunc
            if len(segseqs) > 0 :
                ccs, _ = poa(segseqs, 2, False, 10, -4, -8, -2, -24, -1)
            else :
                ccs = '-'
            eccDNAseq[junction] = ccs
    return revisedJunc, eccDNAseq


def multi_process(data_lst, func, nth, **kargs):
    '''
    multiple processing to handle data list
    :param data_lst: data list
    :param func: unified approach to the processing of various processes
    :param nth: processor number
    :return: tag_infos [list] returned results for all processors
      
    '''
    if nth > 1:
        sub_tasks = split_bins(data_lst, nth)
        pools = __import__('multiprocessing').Pool(nth)
        _func = __import__('functools').partial(func, **kargs)
        tag_infos = pools.map(_func, sub_tasks) # multiple processing
        pools.close(); pools.join()
        #tag_infos = [subline for line in tag_infos for subline in line]
    else:
        tag_infos = func(data_lst, **kargs)
    
#    import pdb; pdb.set_trace()
    #JUNC, DNSeq = [idx[0] for idx in tag_infos ], [idx[1] for idx in tag_infos ]
    return tag_infos



def FullSeqs(Junctions, SAInfo, input_fq, LinearReads, Tag, mapq_cutoff, merge_dist,verbose, threads, filtering_level):
    if verbose >= 3 :
        print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Full sequence\n")
    revisedJunc = {}
    eccDNAseq = {}
    readseq = {}
    for fq in SeqIO.parse(input_fq, "fastq"):
        readseq[fq.name] = str(fq.seq)
    #num = 0
    #total = len(Junctions)
    #if verbose >= 3:
    #    pbar = ProgressBar().start()
    kargs = defaultdict()
    kargs['Junctions'], kargs['Tag'], kargs['SAInfo'], kargs['LinearReads'], kargs['mapq_cutoff'], kargs['merge_dist'], kargs['verbose'] = Junctions, Tag, SAInfo, LinearReads, mapq_cutoff, merge_dist, verbose
    kargs['readseq'], kargs['thread'], kargs['filtering_level'] = readseq, threads, filtering_level
    results = multi_process(
            Junctions.keys()        , 
            consensusSeq   ,
            kargs['thread']     , 
            kargs = kargs
        )
    for thread in range(threads) :
        revisedJunc.update(results[thread][0])
        eccDNAseq.update(results[thread][1])
#    import pdb; pdb.set_trace()
    return(eccDNAseq, revisedJunc)

#eccDNAseq, revisedJunc = FullSeqs(Junctions, SAInfo, input_fq, LinearReads, Tag, mapq_cutoff, merge_dist,verbose, threads)





def Ori_Coverage_count(ont_bam, Junctions, revisedJunc, verbose, filtering_level) :
    '''
    :param ont_bam:
    :param Junctions:
    :return:
    '''
    if verbose >= 3 :
        print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Coverage count around junctions\n")
    if filtering_level == 'high' :
        supporting_reads = 5
    elif filtering_level == 'low' :
        supporting_reads = 2
    else :
        supporting_reads = 1
    bamFile = ps.AlignmentFile("%s" % ont_bam, "rb")
    reference_contigs = bamFile.header['SQ']
    header_dict = {}
    for reference in reference_contigs:
        header_dict[reference['SN']] = reference['LN']
    ratiocovL = {}
    ratiocovR = {}
    mannwhitneyuL = {}
    mannwhitneyuR = {}
    cov_range = 50
    num = 0
    total = len(Junctions)
    if verbose >= 3 :
        pbar = ProgressBar().start()
    if 'chrM' in header_dict:
        covchrm = bamFile.count_coverage(contig='chrM')
        sumchrM = list(np.uint32(np.array([covchrm[0], covchrm[1], covchrm[2], covchrm[3]]).sum(axis=0)))
    for junction in Junctions :
        if len(Junctions[junction]) < supporting_reads :
            num += 1
            if verbose >= 3:
                pbar.update(int((num / (total)) * 100))
            continue
        chrom = junction.split('\t')[0]
        junL = int(revisedJunc[junction].split('\t')[1])
        junR = int(revisedJunc[junction].split('\t')[2])
        if chrom == 'chrM' :
            sumLout = sumchrM[max(junL - cov_range, 0):junL]
            sumLout = [0] * (cov_range - len(sumLout)) + sumLout
            sumLin = sumchrM[junL:min(junL+cov_range, header_dict[chrom])]
            sumLin = sumLin + [0] * (cov_range - len(sumLin))
            sumRout = sumchrM[junR:min(junR+cov_range, header_dict[chrom])]
            sumRout = sumRout + [0] * (cov_range - len(sumRout))
            sumRin = sumchrM[max(junR-cov_range, 0):junR]
            sumRin = [0] * (cov_range - len(sumRin)) + sumRin
        else :
            if junL == 0 :
                sumLout = [0]*cov_range
            else :
                covLout = bamFile.count_coverage(contig=chrom, start=max(junL - cov_range, 0), stop=junL)
                sumLout = list(np.uint32(np.array([covLout[0], covLout[1], covLout[2], covLout[3]]).sum(axis=0)))
                sumLout = [0]*(cov_range-len(sumLout)) + sumLout
            covLin = bamFile.count_coverage(contig=chrom, start=junL, stop=min(junL+cov_range, header_dict[chrom]))
            sumLin = list(np.uint32(np.array([covLin[0], covLin[1], covLin[2], covLin[3]]).sum(axis=0)))
            sumLin = sumLin + [0]*(cov_range-len(sumLin))
            if junR == header_dict[chrom] :
                sumRout = [0]*cov_range
            else :
                covRout = bamFile.count_coverage(contig=chrom, start=junR, stop=min(junR+cov_range, header_dict[chrom]))
                sumRout = list(np.uint32(np.array([covRout[0], covRout[1], covRout[2], covRout[3]]).sum(axis=0)))
                sumRout = sumRout + [0] * (cov_range - len(sumRout))
            covRin = bamFile.count_coverage(contig=chrom, start=max(junR-cov_range, 0), stop=junR)
            sumRin = list(np.uint32(np.array([covRin[0], covRin[1], covRin[2], covRin[3]]).sum(axis=0)))
            sumRin = [0]*(cov_range-len(sumRin)) + sumRin
        ratioL = sum(sumLin) / (sum(sumLout) + sum(sumLin) + 1)
        ratiocovL[junction] = ratioL
        if sumLout == sumLin:
            mannwhitneyuL[junction] = 1
        else:
            mannwhitneyuL[junction] = stats.mannwhitneyu(sumLin, sumLout, alternative='greater').pvalue
        ratioR = sum(sumRin) / (sum(sumRout) + sum(sumRin) + 1)
        ratiocovR[junction] = ratioR
        if sumRout == sumRin:
            mannwhitneyuR[junction] = 1
        else:
            mannwhitneyuR[junction] = stats.mannwhitneyu(sumRin, sumRout, alternative='greater').pvalue
        num += 1
        if verbose >= 3:
            pbar.update(int((num / (total)) * 100))
    return(ratiocovL, ratiocovR, mannwhitneyuL, mannwhitneyuR)

#ratiocovL, ratiocovR, mannwhitneyuL, mannwhitneyuR = Coverage_count(ont_bam, Junctions, revisedJunc, verbose, filtering_level)

def Coverage_count(covfile, Junctions, revisedJunc, verbose, filtering_level) :
    if verbose >= 3 :
        print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Coverage count around junctions\n")
    if filtering_level == 'high' :
        supporting_reads = 5
    elif filtering_level == 'low' :
        supporting_reads = 2
    else :
        supporting_reads = 1
    ratiocovL = {}
    ratiocovR = {}
    mannwhitneyuL = {}
    mannwhitneyuR = {}
    cov_range = 50
    num = 0
    total = len(Junctions)
    if verbose >= 3 :
        pbar = ProgressBar().start()
    Covdict = defaultdict(int)
    with open(covfile,"r") as covFile:
        for line in covFile:
            chrom, pos, depth= line.strip().split("\t")
            Pos = chrom + '\t' + pos
            Covdict[Pos] = int(depth)
    for junction in Junctions :
        if len(Junctions[junction]) < supporting_reads :
            num += 1
            if verbose >= 3:
                pbar.update(int((num / (total)) * 100))
            continue
        chrom = junction.split('\t')[0]
        junL = int(revisedJunc[junction].split('\t')[1])
        junR = int(revisedJunc[junction].split('\t')[2])
        sumLin = []
        sumLout = []
        sumRin = []
        sumRout = []
        for dist in range(cov_range) :
            LoutPos = chrom + '\t' + str(junL - (cov_range - 1 - dist) )
            LinPos = chrom + '\t' + str(junL + dist + 1 )
            RinPos = chrom + '\t' + str(junR - (cov_range - 1 - dist) )
            RoutPos = chrom + '\t' + str(junR + dist + 1  )
            sumLin.append(Covdict[LinPos])
            sumLout.append(Covdict[LoutPos])
            sumRin.append(Covdict[RinPos])
            sumRout.append(Covdict[RoutPos])
        ratioL = sum(sumLin) / (sum(sumLout) + sum(sumLin) + 1)
        ratiocovL[junction] = ratioL
        if sumLout == sumLin:
            mannwhitneyuL[junction] = 1
        else:
            mannwhitneyuL[junction] = stats.mannwhitneyu(sumLin, sumLout, alternative='greater').pvalue
        ratioR = sum(sumRin) / (sum(sumRout) + sum(sumRin) + 1)
        ratiocovR[junction] = ratioR
        if sumRout == sumRin:
            mannwhitneyuR[junction] = 1
        else:
            mannwhitneyuR[junction] = stats.mannwhitneyu(sumRin, sumRout, alternative='greater').pvalue
        num += 1
        if verbose >= 3:
            pbar.update(int((num / (total)) * 100))
    return(ratiocovL, ratiocovR, mannwhitneyuL, mannwhitneyuR)



def output(OnesegJunction, OnesegFa, Junctions, Tag, ratiocovL, ratiocovR, mannwhitneyuL, mannwhitneyuR, Label, filtering_level, eccDNAseq, revisedJunc, verbose) :
    eccDNAnum = 0
    FulleccDNAnum = 0
    if filtering_level == 'high' :
        supporting_reads = 5
        P_value = 0.01
        ratio = 0.6
    elif filtering_level == 'low' :
        supporting_reads = 2
        P_value = 0.05
        ratio = 0.5
    else :
        supporting_reads = 1
        P_value = 1
        ratio = 0
    if verbose >= 3 :
        print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"),
              "Write Single Segments junctions to OutFile\n")
    BreakJunc = []
    FullJunc = []
    outFile = open(OnesegJunction, "w+")
    outFa = open(OnesegFa, "w+")
    for junction in Junctions :
        FullJuncNum = 0
        FullBSJNum = 0
        BreakBSJNum = 0
        for read in Junctions[junction] :
            if Tag[read] == 'FullJunc' :
                FullJuncNum += 1
            elif Tag[read] == 'FullBSJ' :
                FullBSJNum += 1
            elif Tag[read] == 'BreakBSJ' :
                BreakBSJNum += 1
        if FullBSJNum == 0 and FullJuncNum == 0 :
            tag = Label + '\tBreak'
            BreakJunc.append(junction)
        else :
            tag = Label + '\tFull'
            FullJunc.append(junction)
        readscount = FullJuncNum + FullBSJNum + BreakBSJNum
        if (readscount >= supporting_reads) and (ratiocovL[junction] > ratio) and (ratiocovR[junction] > ratio) and (mannwhitneyuL[junction] < P_value) and (mannwhitneyuR[junction] < P_value) :
            juncrecord = '\t'.join([revisedJunc[junction], tag, str(FullJuncNum + FullBSJNum), str(BreakBSJNum), str(mannwhitneyuL[junction]), str(mannwhitneyuR[junction]),
                            str(ratiocovL[junction]), str(ratiocovR[junction]), str(Junctions[junction])])
            outFile.write(juncrecord + '\n')
            eccDNAnum += 1
            if tag == Label + '\tFull' :
                junc = JuncInfo(revisedJunc[junction].split('\t'))
                farecord = '> ' + junc.chrom + ':' + str(junc.startpos) + '-' + str(junc.endpos) + '\t' + tag + '\n' + eccDNAseq[junction] + '\n'
                outFa.write(farecord)
                FulleccDNAnum += 1
    outFa.close()
    outFile.close()
    return(FullJunc, BreakJunc, eccDNAnum, FulleccDNAnum)

#FullJunc, BreakJunc = output(OnesegJunction, OnesegFa, Junctions, Tag, ratiocovL, ratiocovR, mannwhitneyuL, mannwhitneyuR, label, filter_level,eccDNAseq, revisedJunc)




def multisegs_cluster(data, maxgap) :
    data.sort()
    groups = [[data[0]]]
    segnum = len(data[0])
    for x in data[1:]:
        newGroup = True
        for group in groups[1:] :
            inGroup = True
            for seg in range(segnum) :
                if (x[seg][0] == group[-1][seg][0]) and (abs(x[seg][1] - group[-1][seg][1]) <= maxgap) and \
                        (abs(x[seg][2] - group[-1][seg][2]) <= maxgap) :
                    continue
                else :
                    inGroup = False
                    break
            if inGroup == True :
                newGroup = False
                group.append(x)
        if newGroup == True :
            groups.append([x])
    return groups



def Multiple_Segment_Merge(multiseg, merge_dist) :
    print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Merge and Count Multiple Segments junctions\n")
    if len(multiseg) == 0 :
        return True
    else :
        Seginfo = {}
        Segreads = defaultdict(list)
        for segnum in multiseg :
            if len(multiseg[segnum]) < 2 :
                continue
            else :
                allsegs = []
                segreads = defaultdict(list)
                for segments in multiseg[segnum] :
                    minseg = min(segments[0:segnum])
                    segqueue = Queue(segnum)
                    for minnum in range(0,segnum) :
                        if segments[minnum] == minseg :
                            break
                        else :
                            segqueue.put(segments[minnum])
                    sortedSegs = segments[minnum:segnum]
                    for item in range(0,segqueue.qsize()) :
                        sortedSegs.append(segqueue.get())
                    for segment in sortedSegs :
                        segment[1] = int(segment[1])
                        segment[2] = int(segment[2])
                    segreads[str(sortedSegs)].append(segments[segnum])
                    allsegs.append(sortedSegs)
                segsgroup = multisegs_cluster(allsegs, merge_dist)
                for group in segsgroup :
                    Seginfo[str(group[0])] = group
                    for item in group :
                        Key = str(group[0])
                        Segreads[Key].append(segreads[str(item)])
        for junction in Segreads :
            reads = str(Segreads[junction]).replace('[','').replace(']','').replace("'",'').split(', ')
            Segreads[junction] = list(np.unique(reads))
    return(Seginfo, Segreads)





def MS_PseudoReference(OutDir, Segreads, refFa, fastq, threadnum, verbose) :
    if verbose >= 3 :
        print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"),
              "Generate PseudoReference for Multiple-segments eccDNA\n")
    ### mkdir temp file
    if not os.path.exists(OutDir + '/MStemp'):
        os.makedirs(OutDir + '/MStemp')
    tmpbed = OutDir + '/MStemp/' + "temp.multiseg.bed"  ## bed file for all segments of mulitiple-fragments eccDNA
    tmpfa = OutDir + '/MStemp/' + "temp.multiseg.fa"  ## fasta file for all segments of mulitiple-fragments eccDNA
    tmpid = OutDir + '/MStemp/' + "temp.multiseg.read.id"  ## full-length reads of mulitiple-fragments eccDNA
    tmpPRfa = OutDir + '/MStemp/' + "temp.PseudoReference.fa"  ## PseudoReference of mulitiple-fragments eccDNA
    tmpfq = OutDir + '/MStemp/' + "temp.multiseg.read.fastq"
    tmpsam = OutDir + '/MStemp/' + "temp.2PseudoReference.sam"
    tmpbam = OutDir + '/MStemp/' + "temp.2PseudoReference.sorted.bam"
    tmpbai = OutDir + '/MStemp/' + "temp.2PseudoReference.sorted.bam.bai"
    with open(tmpbed, "w+") as TEMPBED:
        with open(tmpid, "w+") as TEMPID:
            for junction in Segreads:
                TEMPBED.write(
                    junction.replace("'", "").replace("], [", "\n").replace(", ", "\t").replace("[", "").replace("]","") + '\n')
                for read in Segreads[junction]:
                    TEMPID.write(read + '\n')
    get_segfa = subprocess.call(["bedtools", "getfasta", "-fi", refFa, "-bed", tmpbed, "-fo", tmpfa], shell=False)  ##
    if get_segfa != 0:
        print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"),
              "An error happened during bedtools getfasta for PseudoReference generation. (bedtools can't read fasta index in presence of hla alts.) Exiting")
        #sys.exit()
    seqs = {}
    with open(tmpfa, 'r') as fa:
        for record in SeqIO.parse(fa, 'fasta'):
            seqs[record.id] = record.seq
    ecclist = []
    for record in Segreads:
        eccdna = record.strip().lstrip('[[').rstrip(']]').split('], [')
        eccid = ''
        pseudoref = ''
        for segment in eccdna:
            segment = segment.split(', ')
            segid = segment[0].strip("'") + ':' + str(segment[1]) + '-' + str(segment[2])
            eccid += segid + '|'
            pseudoref += seqs[segid]
        eccid = eccid.rstrip('|')
        eccref = SeqRecord(pseudoref, id=eccid)
        ecclist.append(eccref)
    SeqIO.write(ecclist, tmpPRfa, "fasta")
    with open(tmpfq, "wb") as TMPFQ:
        get_multireads = subprocess.call(
            ["seqtk", "subseq", fastq, tmpid], stdout=TMPFQ)
    if get_multireads != 0:
        print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"),
              "An error happened during seqtk subseq for full-length reads extraction. Exiting")
        #sys.exit()
    realignment = subprocess.call(["minimap2", "-t", str(threadnum), "-ax", "map-ont", tmpPRfa, tmpfq, "-o", tmpsam],
                                  shell=False)
    samtoolsSort = subprocess.call(["samtools", "sort", "-@", str(threadnum), "-O", "bam", "-o", tmpbam, tmpsam],
                                   shell=False)
    samtoolsIndex = subprocess.call(["samtools", "index", "-@", str(threadnum), tmpbam, tmpbai], shell=False)
    removetemp = subprocess.call(["rm", "-f", tmpbed, tmpfa, tmpid, tmpPRfa, tmpsam], shell=False)
    if (realignment + samtoolsSort + samtoolsIndex) != 0:
        print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"),
              "An error happened during minimap2 (or samtools) for full-length reads alignment. Exiting")
        #sys.exit()
    return (tmpbam, tmpfq)
#tmpbam, tmpfq = MS_PseudoReference(OutDir, Segreads, refFa, fastq, threadnum, verbose)
#multibegin = time.time()
#MS_SAInfo, MSLinearReads, MSreadlens, MSreadnum, MSmappedreadsnum = get_continuous_SAinfo(tmpbam,read_gap, mapq_cutoff, verbose, multibegin)


def parse_junction_from_SAinfo_of_MultipleSegment_FLreads(SAInfo, readlens, verbose) :
    if verbose >= 3:
        print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Parse continuous SAinfo  for multiple-segments eccDNA detection\n")
    Tag = {}  # Full, Break
    Juncinfo = {}
    Junctions = defaultdict(list)
    for readid in SAInfo.keys():
        if len(SAInfo[readid]) < 2 :
            Tag[readid] = "Break"
            #print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"),"An error happenend during execution. Exiting")
            #sys.exit()
        else :
            info = SAInfo[readid]
            SAnum = len(info)
            DG = nx.DiGraph()
            lastnode = info[0][2]
            for subscript in range(1, SAnum):
                newnode = info[subscript][2]
                if DG.has_edge(lastnode, newnode):
                    DG.edges[lastnode, newnode]['weight'] += 1
                else:
                    DG.add_edge(lastnode, newnode, weight=1)
            cycles = list(nx.simple_cycles(DG))
            maxweight = 0
            for cycle in cycles:
                if len(cycle) == 1:
                    Tag[readid] = "Full"
                    cycleweight = DG.get_edge_data(cycle[0], cycle[0])['weight']
                    if cycleweight > maxweight:
                        maxweight = cycleweight
                        bestcycle = cycle[0]
            if maxweight == 0 :
                Tag[readid] = 'invalid'
            else:
                cyclelen = 0
                for item in info :
                    if item[2] == bestcycle :
                        cyclelen += (item[1]-item[0])
                if cyclelen / readlens[readid] >= 0.5:
                    Juncinfo[readid] = bestcycle
                    Junctions[bestcycle].append(readid)
                    Tag[readid] = 'Full'
                else :
                    Tag[readid] = 'invalid'
    return (Juncinfo, Junctions, Tag)
#MS_Juncinfo, MS_Junctions, MS_Tag = parse_junction_from_SAinfo_of_MultipleSegment_FLreads(MS_SAInfo, MSreadlens, verbose)


def FullSeqs_MS(MS_Junctions, MS_SAInfo, input_fq, verbose, filtering_level):
    if verbose >= 3 :
        print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"), "Full sequences for multiple-segments eccDNA\n")
    eccDNAseq = defaultdict(str)
    if filtering_level == 'high' :
        supporting_reads = 5
    elif filtering_level == 'low' :
        supporting_reads = 2
    else :
        supporting_reads = 1
    readseq = {}
    for fq in SeqIO.parse(input_fq, "fastq"):
        readseq[fq.name] = str(fq.seq)
    for junction in MS_Junctions.keys():
        if len(MS_Junctions[junction]) < supporting_reads :
            ccs = "-"
        else:
            segseqs = []
            CCS = False
            for readid in MS_Junctions[junction] :
                if len(MS_SAInfo[readid]) > 2 :
                    CCS = True
                    sainfo = MS_SAInfo[readid]
                    for sanum in range(len(sainfo) - 1):
                        segment = AlignSeg(sainfo[sanum])
                        nextseg = AlignSeg(sainfo[sanum + 1])
                        if (junction == segment.chrom) and (junction == nextseg.chrom) :
                            misbase = nextseg.readS - segment.readE
                            seq = readseq[readid][segment.readS:nextseg.readS]
                            segseqs.append(seq)
                    sanum += 1
                    segment = AlignSeg(sainfo[sanum])
                    if (junction == segment.chrom) :
                        seq = readseq[readid][segment.readS:segment.readE]
                        segseqs.append(seq)
            #print(junction + '\t' + str(len(segseqs)))
            if (not CCS) :
                ccs = "-"
            else:
                ccs, _ = poa(segseqs, 2, False, 10, -4, -8, -2, -24, -1)
            #print(junction + '\t' + str(len(segseqs)) + '\t' + str(len(ccs)))
        eccDNAseq[junction] = ccs
    return (eccDNAseq)
#MSeccDNAseq = FullSeqs_MS(MS_Junctions, MS_SAInfo, tmpfq, verbose, filtering_level)



def MSoutput(MulsegJunction, MulsegFa, MS_Junctions, MSeccDNAseq, Label, filtering_level, verbose) :
    FullMSeccDNAnum = 0
    if filtering_level == 'high':
        supporting_reads = 5
        P_value = 0.01
        ratio = 0.6
    elif filtering_level == 'low':
        supporting_reads = 2
        P_value = 0.05
        ratio = 0.5
    else:
        supporting_reads = 1
        P_value = 1
        ratio = 0
    if verbose >= 3:
        print(datetime.datetime.now().strftime("\n%Y-%m-%d %H:%M:%S:"),
              "Write Multiple Segments junctions to OutFile\n")
    with open(MulsegJunction, "w+") as outFileMS :
        with open(MulsegFa, "w+") as outFaMS:
            for junction in MS_Junctions.keys():
                if len(MS_Junctions[junction]) >= supporting_reads :
                    segments = junction.split("|")
                    juncrecord = '\t'.join([junction, Label, "Full", str(len(segments)), str(len(MS_Junctions[junction]))])
                    outFileMS.write(juncrecord + '\n')
                    farecord = '> ' + junction + '\t' + Label + '\tFull' + '\n' + MSeccDNAseq[junction] + '\n'
                    outFaMS.write(farecord)
                    FullMSeccDNAnum  += 1
    return (FullMSeccDNAnum)
#FullMSeccDNAnum = MSoutput("MS_OUT.out", "MS_FA.fa", MS_Junctions, MSeccDNAseq, Label, filtering_level, verbose)




