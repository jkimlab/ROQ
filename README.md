# ROQ
Next-generation sequencing (NGS) has advanced genomic research, highlighting the need for accurate read alignment metrics. The Read Overlapping Quality (ROQ) score is a novel machine-learning-based metric that predicts the overlap between a read’s true origin and its mapped position. Unlike conventional metrics, ROQ better reflects true read origins, even in complex genomic regions, enhancing downstream analyses by reducing misalignment and improving read quality assessment.


## Generate ROQ with example data
```bash
# Clone the repository
git clone https://github.com/jkimlab/ROQ/
cd ./ROQ/example/

#Feature extraction
perl ../bin/ExtractFeatures.pl example.sam  1> feature.txt 2> log

#Predict ROQ
python ../bin/predROQ.py -i feature.txt -ib example.sam -m ../model/hg38_bowtie_model.pkl -o ./ 

# View the updated BAM file
samtools view output.bam  

```

## Commands
```bash
#Feature extraction
perl [path to ROQ bin]/ExtractFeatures.pl [user input read alignment]  1> feature.txt 2> log

#Predict ROQ
python [path to ROQ bin]/predROQ.py -i feature.txt -ib [user input read alignment] -m [path to ROQ model]/hg38_bowtie_model.pkl -o [output directory] 

```


## Input - Output example
- Input: `.SAM` or `.BAM` formatted read alignment
- `feature.txt`: extracted features for each read
  ```text
  READ_NAME	MAPQ	ALIGN_SCORE	SECONDARY_ALIGN_SCORE	MISMATCHES	GAP_OPENS	GAP_EXT	EDIT_DIST	MATE_ALIGNMENT_SCORE	ALIGN_SCORE_DIFF	INSERT_SIZE	READ_GC_CONT	N_LOW_QUALITY_BASE	AVG_QUALITY_BASE_SCORE
  1-26020508	30	0	-5	0	0	0	0	0	5	200	0.533	0	69.540
  1-26020508	30	0	-5	0	0	0	0	0	5	-200	0.513	0	70.140
  1-26020324	33	-2	-90	1	0	0	1	0	88	192	0.407	0	68.387
  1-26020324	33	0	-78	0	0	0	0	-2	78	-192	0.407	0	69.607
  1-26020524	27	-2	-80	1	0	0	1	0	78	193	0.373	0	69.500
  1-26020524	27	0	-73	0	0	0	0	-2	73	-193	0.420	0	69.280
  ```
- Output: `output.bam`, ROQ tag appended read alignment in `.BAM` format file.
  ```text
  #results of
  $samtoole view output.bam
  
  1-26020332      99      1       108376487       42      150M    =       108376538       201     TTTTATTGTGTACAGGTGTGTTTTCTATTCCTGCTATCTCCATGACTTTTAAGATGAGGGGACATTAGAATGTGTCAAAGGCTTTTTTCAGCATCTAGTGAAATGATCATGCATTAGTGTGTGTGTGTGAGTGTTTGTGTGTTCTGTAAG   =CCGGGGGGGGGGJCJJJGJJJJJJJJGJJJJGJGJJGGJJJJGJJJJGCGJCCG=JCGJJJGCGJGC=JCGGJGGGGCGGGCGGG=GGCGG8=GGGGGGJGGCGCGGGG=GGGGCCG=CCGGCGGGGGGGGC=(=CGGCGGCCC=GGGC   AS:i:-2 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:134G15     YS:i:-3 YT:Z:CP RQ:f:1.29227
  1-26020332      147     1       108376538       42      150M    =       108376487       -201    AGATGAGGGGACATTAGAATGTGTCAAAAGCTTTTTTCAGCATCTAGTGAAATGATCATGCATTAGTGTGTGTGTGTGAGTGTGTGTGTGTTCTGTAAGTTTCTTTATATGGTGGATTAAATTGTACTTTTGTATGTTGAACCATTCCTG   GGGGGC1GGG=G=GGCGGGCGGGGGGCG1GCCG8G1GCGGGCCJ=JJCGGGCGGCGGGCGG=GGGCC8CCCGCC1GCGGGCGCCGG=JCCJJGGGGJCJG=CJJJCJJJCCCGJCJGJJGJGJJJGJJCJJJGJJJJGGGGGGGGGG==C   AS:i:-3 XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:28G121     YS:i:-2 YT:Z:CP RQ:f:1.29227
  ```

## ROQ tag (Read Overlap Quality tag, represented as RQ:f:value)
The ROQ tag, which stands for Read Overlap Quality, is added to each read in the BAM file by this tool. This custom tag represents a quantitative metric that evaluates how accurately a read overlaps with its true genomic origin. It serves as a more reliable alternative to traditional alignment quality metrics like MAPQ.

### Details of the RQ Tag
- **Full Name**: Read Overlap Quality
- **Type**: Floating-point value
- **Purpose**: Quantify the overlap between the read's mapped position and its true origin.
- **Key Benefit**: Provides a robust metric for filtering misaligned reads or assessing alignment quality.
- **Example**: `RQ:f:1.29227`
  



## Requirements

### System Requirements

| Software          | Version       | Installation Guide                                    |
|--------------------|---------------|------------------------------------------------------|
| **Perl**          | 5.30.0 or higher | [Perl Installation Guide](https://www.perl.org/get.html) |
| **Python**        | 3.8.12 or higher | [Python Downloads](https://www.python.org/downloads/)  |
| **Samtools**      | 1.9  | [Samtools GitHub](https://github.com/samtools/samtools) |
| **convert2bed**   | 2.4.38        | [bedops Github](https://github.com/bedops/bedops)     |

### Required Perl Modules
- `Cwd` (Core module)
- `FindBin` (Core module)

### Required Python Libraries
| Library         | Version   | Installation Command               |
|------------------|-----------|------------------------------------|
| `pandas`        | 1.5.3     | `pip install pandas==1.5.3`        |
| `numpy`         | 1.23.0    | `pip install numpy==1.23.0`        |
| `pysam`         | 0.22.0    | `pip install pysam==0.22.0`        |
| `scikit-learn`  | 1.2.2     | `pip install scikit-learn==1.2.2`  |
| `xgboost`       | 1.7.6     | `pip install xgboost==1.7.6`       |


## Contact
ibclab.kr@gmail.com
