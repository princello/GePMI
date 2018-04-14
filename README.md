# GePMI

GePMI: A statistical model for personal intestinal microbiome identification

(Generating inter-individual similarity distribution for Personal Microbiome Identification)

## Description
GePMI is a framework for judging metagenomic samples from intra or inter-individual.

The principle of GePMI is based on pairwise similarities between the metagenomes of any two individuals obey a Beta distribution and that a p-value derived accordingly well characterizes whether two samples are from the same individual or not. To control the false discovery rate (FDR) in multiple testing, Benjamini and Yekutieli’s method was used to transform p–values to q-values. So GePMI can help you to determine whether the two samples are from the same individual through three thresholds: similarities of input, p-values and q-values.


## Quick Tutorial
### Input file 
Similarity matrix，like `test-100-18-10000.csv`

MRA_P1E-0 | MRA_P1E-7| MRA_P1E-90 | ......
------------- | ------------- | ------------- | -------------
1 | 0.2435 | 0.1969 | ......
0.2435 | 1 | 0.2130 | ......
0.1969 | 0.2130 | 1 | ......
...... | ...... | ...... | ......

### Running script
`python GePMI.py -i test-100-18-10000.csv -o output`
### Output files
`result.txt`: Identification of Test sample and Target sample from one individual under threshold

The first column shows the target samples used to generate the distribution;
The following columns are test samples which are significantly similar to the target samples，That is to say, GePMI judgment this(theses) test(s) and target sample come from the same individual

A glimpse of the result.txt:

Target_sample | Significant_similar_test_sample 1| Significant_similar_test_sample 2 | ......
------------- | ------------- | ------------- | -------------
MRA_P1E-0| MRA_P1E-7 | MRA_P1E-90 |
MRA_P1E-7 | MRA_P1E-0 | MRA_P1E-90 |
MRA_P1E-90 | No matches 
...... |

`detail.txt`: Values(MinHash similarity, GePMI p_values, GePMI q_values) of each Test sample in Target sample inter-individual similarity distribution

### A brief manual:
![image](https://github.com/princello/GePMI/blob/master/brief%20manual.jpg)

### Contact
Please report any problems directly to the GitHub issue tracker

Also, you can send your feedback to wangzicheng12@mails.tsinghua.edu.cn

-

## Details
### 0.Preparation(Optional)
down sample to the same size

### 1.Input
[sourmash](https://github.com/dib-lab/sourmash) csv file

#### Recommended usage of sourmash
##### From sequence to signature
`sourmash compute -k 18 -n 10000 subject-sample.fa -o subject-sample.sig`
>for `subject-sample`,if you know the subject name and sample number,for example:
>>HMP_subject1-v1.fa, please name it like this
##### Calculating similarity between signatures
`sourmash compare -k 18 --csv test-100-18-10000.csv -o test *.sig`
>We suggest that you name csv output file like this
>>test-100-18-10000.csv, -18 is KSIZES, -10000 is NUM_HASHES used in sourmash

### 2.GePMI Usage

`python GePMI.py -i input.csv -p 0.001 -q 0.01 -s 0 -o outputDir -t`

-i your input sourmash csv file

>Please name it as follows:
>>[prefix]-[base number(millions)]-[k-mer length]-[hashes used in sourmash].csv

>for eaxamle:
>>test-100-18-10000.csv

* -p threshold of p-value (default 0.001)

* -q threshold of q-value (default  0.01)

* -s threshold of similarity (default  0)

* -o path of output dirctory (default ./output)

* -t save p/q value matrix or not


* Figures (-t for generating)

roc-xx-xx-xx.pdf  : Receiver Operating Characteristic of known sample for PMI

prc-xx-xx-xx.pdf  : Precesion-Recall Curve of known sample for PMI

fdr-xx-xx-xx.pdf  : False Discovery Rate Comparing of p/q values

fdr-xx-xx-xx.d.pdf: Detail of FDR plot in range of 0 - 0.01

Files (-t for generating)

p_values.txt : p value of each test sample in target samples's inter-individual similarity distribution

q_values.txt : q value by Benjamini & Yekutieli method (Benjamini Y, Yekutieli D. The control of the false discovery rate in multiple testing under dependency[J]. Annals of Statistics, 2001, 29(4):1165--1188.)
