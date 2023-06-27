# FLED: a full-length eccDNA detector for long-reads se-quencing data
Full Length eccDNA detection
Reconstructing the full-length sequence of extrachromosomal circular DNA (eccDNA) from short sequencing reads has proved challenging given the similarity of eccDNAs and their corresponding linear DNAs. Previous sequencing methods were unable to achieve high-throughput detection of full-length eccDNAs. Herein, a novel algorithm was developed, called Full-Length eccDNA Detection (FLED), to reconstruct the sequence of eccDNAs based on the strategy that combined rolling circle amplification (RCA) and nanopore long-reads sequencing technology. Seven human epithelial and cancer cell line samples were analyzed by FLED and over 5,000 full-length eccDNAs were identified per sample. The structures of identified eccDNAs were validated by both PCR and Sanger sequencing. Compared to other published nanopore-based eccDNA detectors, FLED exhibited higher sensitivity. In cancer cell lines, the genes overlapped with eccDNA regions were enriched in cancer-related pathways and cis-regulatory elements can be predicted in the upstream or downstream of intact genes on eccDNA molecules, and the expressions of these cancer-related genes were dysregulated in tumor cell lines, indicating the regulatory potency of eccDNAs in biological processes. The proposed method takes advantage of nanopore long reads and enables unbiased reconstruction of full-length eccDNA sequences.

## Requirements
* Software
* Package (python 3 +)
    - [pysam](http://pysam.readthedocs.org/en/latest/) (>=0.16)
    - [networkx](https://github.com/networkx/networkx)(>=2.5)
    - [progressbar](https://pypi.org/project/progressbar)(>=2.5)
    - [biopython](https://biopython.org/)(>=1.76)
    - [numpy](https://numpy.org/)(>=1.19.1)
    - [pyspoa](https://pypi.org/project/pyspoa/)(==0.0.6)
    - [scipy](https://pypi.org/project/scipy/)(==1.5.3)
 
## Installation
```bash
git clone https://github.com/FuyuLi/FLED.git
cd FLED
python ./setup.py install
```

## Required files
Users can prepare the external files under the following instructions:
1) Indexed genome fasta file
```
samtools faidx $genome
```

2) Indexed coordinate sorted bam file
```
minimap2 -t 10 -ax map-ont $genome $fastq | samtools sort -@ 4 -O bam -o align.sorted.bam
samtools index  -@ 10 align.sorted.bam align.sorted.bam.bai
```

## Usage
```
usage: FLED <subprogram> [options]

Commands:
   Detection       Identify extrachromosomal circular DNA from Nanopore reads
   Simulate        Simulate extrachromosomal circular DNA

```


### Detection
```
usage: FLED Detection [options]

Identify extrachromosomal circular DNA from Nanopore reads

required arguments:
  -i                    Input: Indexed coordinate sorted bam file
  -fq                   Input: Nanopore fastq file

optional arguments:
  -o , --outPrefix          Prefix of Output filename
  -dir , --directory        Output directory, default is the working directory
  -q , --mapq               Minimum mapping quality allowed in the alignments. Default: 10
  -g , --gap                Maximum gap allowed in the Continuous supplementary alignments reconstruction. Default: 50
  -K , --clustering_dist    Cluster ecDNAs that are K nucleotides apart in the same interval. Default: 50
  -v , --verbose            Verbose level, 1=error,2=warning, 3=message. Default: 3
  -F , --filter_level       Filtering level, high=5Jreads&HighCoverage, low=2Jreads&MiddleCoverage, none=no-filtering. Default: low
  -t , --threads            Number of threads to use.Default 1
```

### Example
```
FLED Detection -i align.sorted.bam -fq example.q7.fastq -o example -dir FLEDoutput -t 8
```

## Output
| File name         |  Details | 
|   :---            | ---        |
| out.DiGraph.OnesegJunction.fa         | full-length sequences of Full eccDNAs detected by FLED |
| out.DiGraph.OnesegJunction.out        | list of eccDNAs detected |

### out.DiGraph.OnesegJunction.out
| No. | Column name     |  Details | 
|:---:|   :---          | ---        |
|  1  | chrom           | chromosome |
|  2  | start           | start coordinate of eccDNA |
|  3  | end             | end coordinate of eccDNA |
|  4  | strand          | strand of eccDNA |
|  5  | label           | Prefix of Output filename |
|  6  | Tag             | Tag of eccDNA: Full for full-length and Break for breakpoint |
|  7  | Nfullpass       | number of reads contained at least one complete copy of eccDNA |
|  8  | Nbreakpoint     | number of reads covered breakpoint of eccDNA only |
|  9  | L_Pvalue        | Pvalue of Wilcoxon test based on the coverage around Left breakpoint |
|  10 | R_Pvalue        | Pvalue of Wilcoxon test based on the coverage around Right breakpoint |
|  11 | L_covRatio      | coverage ratio inside and outside the Left breakpoint |
|  12 | R_covRatio      | coverage ratio inside and outside the Right breakpoint |
|  13 | readID          | ID of reads supporting this eccDNA |

## Citation
Fuyu L, Wenlong M, Wenxiang L, Ying W, Xiaohan L, Xianjun D, Yunfei B. (2023): FLED: a full-length eccDNA detector for long-reads se-quencing data. bioRxiv. Preprint. 
doi: https://doi.org/10.1101/2023.06.21.545840

## License
