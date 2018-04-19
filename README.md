# GePMI

GePMI: A statistical model for personal intestinal microbiome identification

(Generating inter-individual similarity distribution for Personal Microbiome Identification)

## Description
*  Purpose 
It is well accepted that there is a large number of variants in the human microbiomes, which results in great differences within inter-individuals under the same condition. Here, we aim to identity intra-individual samples through metagenomic similarity.  

When you have a set of metagenomic samples, and you want to know which samples come from the same individual. Or if you want to know whether the intra-individual metagenomic samples are still similar before and after some disturbances, it is a good choice to use GePMI.

The principle of GePMI is based on pairwise similarities between the metagenomes of any two individuals obey a Beta distribution and that a p-value derived accordingly well characterizes whether two samples are from the same individual or not. To control the false discovery rate (FDR) in multiple testing, Benjamini and Yekutieli’s method was used to transform p–values to q-values. So GePMI can help you to determine whether the two samples are from the same individual through three thresholds: similarities of input file, GePMI p-values and GePMI q-values.


## Quick Tutorial
### Input file 
If you want to start analyzing from the initial metagenomic sequence, please jump to details.

The similarity matrix is the initial input of the method. You can generate it in many ways, but we recommend using [sourmash](https://github.com/dib-lab/sourmash), in a word, your input should be arranged as follows:

MRA_P1E-0 | MRA_P1E-7| MRA_P1E-90 | ......
------------- | ------------- | ------------- | -------------
1 | 0.2435 | 0.1969 | ......
0.2435 | 1 | 0.2130 | ......
0.1969 | 0.2130 | 1 | ......
...... | ...... | ...... | ......

The first row is the name of samples, the rest part is the corresponding similarities.

### Running script
`python GePMI.py -i test-100-18-10000.csv -o output`

GePMI.py is the main python script for judging metagenomic samples from intra or inter-individual.

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

`detail.txt`: Values(MinHash similarity, GePMI p-values, GePMI q-values) of each Test sample in Target sample inter-individual similarity distribution

### A brief manual:
![image](https://github.com/princello/GePMI/blob/master/brief%20manual.jpg)

### Contact
Please report any problems directly to the GitHub [issue tracker](https://github.com/princello/GePMI/issues)

Also, you can send your feedback to wangzicheng12@mails.tsinghua.edu.cn


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

* -i your input sourmash csv file

>Please name it as follows:
>>[prefix]-[base number(millions)]-[k-mer length]-[hashes used in sourmash].csv

>for eaxamle:
>>test-100-18-10000.csv

* -p threshold of p-value (default 0.001)

* -q threshold of q-value (default  0.01)

* -s threshold of similarity (default  0)

* -o path of output dirctory (default ./output)

* -t save p/q value matrix or not
>p_values.txt : p value of each test sample in target samples's inter-individual similarity distribution

>q_values.txt : q value by Benjamini & Yekutieli method (Benjamini Y, Yekutieli D. The control of the false discovery rate in multiple testing under dependency[J]. Annals of Statistics, 2001, 29(4):1165--1188.)

* Figures (-t for generating)

>roc-xx-xx-xx.pdf  : Receiver Operating Characteristic of known sample for PMI

>prc-xx-xx-xx.pdf  : Precesion-Recall Curve of known sample for PMI

>fdr-xx-xx-xx.pdf  : False Discovery Rate Comparing of p/q values

>fdr-xx-xx-xx.d.pdf: Detail of FDR plot in range of 0 - 0.01



