# FLED
Full Length eccDNA detection

## Installation requirements
* Software
* Package (python 3 +)
    - [pysam](http://pysam.readthedocs.org/en/latest/) (>=0.16)
    - [networkx](https://github.com/networkx/networkx)(>=2.5)
    - [progressbar](https://pypi.org/project/progressbar)(>=2.5)
    - [biopython](https://biopython.org/)(>=1.76)
    - [numpy](https://numpy.org/)(>=1.19.1)
    - [pyspoa](https://pypi.org/project/pyspoa/)(==0.0.5)
    - [scipy](https://pypi.org/project/scipy/)(>=1.5.3)
 
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
minimap2 -t 10 -ax map-ont $$genome $fastq > align.sam
samtools sort -@ 4 -O bam -o align.sorted.bam align.sam
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

## License
