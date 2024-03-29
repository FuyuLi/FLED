# FLED: a full-length eccDNA detector for long-reads sequencing data

Reconstructing the full-length sequence of extrachromosomal circular DNA (eccDNA) from short sequencing reads has proved challenging given the similarity of eccDNAs and their corresponding linear DNAs. Previous sequencing methods were unable to achieve high-throughput detection of full-length eccDNAs. Herein, a novel algorithm was developed, called Full-Length eccDNA Detection (FLED), to reconstruct the sequence of eccDNAs based on the strategy that combined rolling circle amplification (RCA) and nanopore long-reads sequencing technology. Seven human epithelial and cancer cell line samples were analyzed by FLED and over 5,000 full-length eccDNAs were identified per sample. The structures of identified eccDNAs were validated by both PCR and Sanger sequencing. Compared to other published nanopore-based eccDNA detectors, FLED exhibited higher sensitivity. In cancer cell lines, the genes overlapped with eccDNA regions were enriched in cancer-related pathways and cis-regulatory elements can be predicted in the upstream or downstream of intact genes on eccDNA molecules, and the expressions of these cancer-related genes were dysregulated in tumor cell lines, indicating the regulatory potency of eccDNAs in biological processes. The proposed method takes advantage of nanopore long reads and enables unbiased reconstruction of full-length eccDNA sequences.

## Requirements
* Software
    - [minimap2](https://github.com/lh3/minimap2)
    - [samtools](https://github.com/samtools/samtools)(==1.10 )
    - [bedtools](https://bedtools.readthedocs.io/en/latest/index.html)
    - [seqtk](https://github.com/lh3/seqtk)
* Package (python 3 +)
    - [pysam](http://pysam.readthedocs.org/en/latest/) (==0.22)
    - [networkx](https://github.com/networkx/networkx)(>=2.5)
    - [progressbar](https://pypi.org/project/progressbar)(>=2.5)
    - [biopython](https://biopython.org/)(==1.76)
    - [numpy](https://numpy.org/)(>=1.19.1)
    - [pyspoa](https://pypi.org/project/pyspoa/)(==0.0.6)
    - [scipy](https://pypi.org/project/scipy/)(==1.5.3)
    - [tqdm](https://pypi.org/project/tqdm/)(>=4.51.0)
 
## Installation
```bash
git clone https://github.com/FuyuLi/FLED.git
cd FLED
python ./setup.py install
```
or
```bash
git clone https://github.com/FuyuLi/FLED.git
cd FLED/dist
pip install FLED-1.7.0.tar.gz
```

## Required files
Users can prepare the external files under the following instructions:
1) Indexed genome fasta file
```
samtools faidx $genome
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
  -ref                  Input: reference genome fasta file
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
FLED Detection -ref genome.fa -fq example.q7.fastq -o example -dir FLEDoutput -t 8
```

## Output
| File name         |  Details | 
|   :---            | ---        |
| out.DiGraph.MulsegFullJunction.fa         | full-length sequences of Full multifragment eccDNAs detected by FLED |
| out.DiGraph.MulsegFullJunction.out        | list of full-length multifragment eccDNAs detected |
| out.DiGraph.OnesegJunction.fa         | full-length sequences of Full simple eccDNAs detected by FLED |
| out.DiGraph.OnesegJunction.out        | list of simple eccDNAs detected |

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
Li F, Ming W, Lu W, Wang Y, Li X, Dong X, Bai Y. FLED: a full-length eccDNA detector for long-reads sequencing data. Brief Bioinform. 2023 Sep 22;24(6):bbad388. doi: 10.1093/bib/bbad388. PMID: 37930031; PMCID: PMC10632013.

## License
