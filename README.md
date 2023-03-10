# Constructing the longest contig containing query sequence from given sequencing reads

## Aim
To construct the longest possible cogtig from given reads containing the query sequence
## Description
The program is divided into 2 components:
1. Sequence alignment. Apply Basic Local Alignment Tool (BLAST) preprocessing and extend matched seeds using Smith-Waterman algorithm to perform local alignment, then subset reads with good alignment quality (seeds found). 
2. Sequence assembly. Construct a De Bruijn Graph from filtered reads from the alignment step. Then traceback the graph to form the query sequence by following Eulerian paths. Finally, build out from the assembled sequence to get the final longest sequence.

## Installation
- [python >= 3.8.3](https://www.python.org/downloads/)
- [numpy](https://numpy.org/install/)
- [itertools](https://docs.python.org/3/library/itertools.html)
- [re](https://docs.python.org/3/library/re.html)
- [sys](https://docs.python.org/3/library/sys.html) 

## Command

To run this program in the command line interface type: 
```text
python3 BLASTpreprocess_run.py -QUERYfasta "data/toy_query.fasta" --READSfasta "toy_READS.fasta" --kmer 3 --outputFile "output/stats_output.txt"
```
### Arguments
```text
usage: BLASTpreprocess_run.py [-h] [--QUERYfasta QUERYfasta] [--READSfasta READSfasta] [--kmer kmer]
                              [--outputFile outputFile]

Reconstruct sequence from reads and query for substring PartI: BLAST preprocessing and exact local alignment

optional arguments:
  -h, --help            show this help message and exit
  --QUERYfasta QUERYfasta
                        input file containing query sequene
  --READSfasta READSfasta
                        input file containing database sequences i.e. READS
  --kmer kmer           length of node string
  --outputFile outputFile
                        name of file to write intermediatestatistics output to
```
### Inputs
Query sequence FASTA file example `toy_query.fasta`:
```text
>queryseq
TGTTACGG
```

Reads FASTA file example `toy_READS.fasta`:
```text
>dbseq1
GGTTGACTA
>dbseq2
ATGGATTGCACGCAG
>dbseq3
CGCAGGTTCTCCGGCCGCTT
>dbseq4
CGCTTGGGTGGAG
>dbseq5
AAGTAAAAAAAAAATAAAGC
```

### Output (Part I: BLAST preprocessing and local alignment)
Example intermediate output from part I `stats_output.txt`: 
```text
Query sequence: TGTTACGG

Exploratory alignment results: 
-----------------------------
dbseq1: GGTTGACTA
Length of database sequence: 9
Time spent searching for seeds: 0.00017118453979492188 seconds
kmer and numeric key: GTT,47
Seed position: 3,2
Time Smith-Waterman: 0.0002923011779785156 seconds
percent matched: 44.44
Smith-Waterman max score: 4
G-T-TAC
G T  AC
GTTG-AC
-----------------------------
dbseq2: ATGGATTGCACGCAG
Length of database sequence: 15
Time spent searching for seeds: 0.00021791458129882812 seconds
kmer and numeric key: ACG,6
Seed position: 0,10
Time Smith-Waterman: 0.0004799365997314453 seconds
percent matched: 44.44
Smith-Waterman max score: 3
G-T-TACG
  T  ACG
-TTGCACG
-----------------------------
...
```

### Outputs (Part II: Sequence assembly from the potential matches)
** To be updated **
