# GePMI

GePMI: A statistical model for personal intestinal microbiome identification

GePMI (Generating inter-individual similarity distribution for Personal Microbiome Identification)

## Input
[sourmash](https://github.com/dib-lab/sourmash) csv file

### Recommended usage of sourmash
`sourmash compute -k 18 -n 10000 subject-sample.fa -o subject-sample.sig`
>for `subject-sample`,if you know the subject name and sample number,for example:
>>HMP_subject1-v1.fa, please name it like this

`sourmash compare -k 18 --csv test-100-18-10000.csv -o test *.sig`
>We suggest that you name csv output file like this
>>test-100-18-10000.csv, -18 is KSIZES, -10000 is NUM_HASHES used in sourmash

###



## Usage

`python GePMI.py -i input.csv -p 0.001 -q 0.01 -s 0 -o outputDir -t`

-i your input sourmash csv file

>Please name it as follows:
>>[prefix]-[base number(millions)]-[k-mer length]-[hashes used in sourmash].csv

>for eaxamle:
>>test-100-18-10.csv

-p threshold of p-value (default 0.001)

-q threshold of q-value (default  0.01)

-s threshold of similarity (default  0)

-o path of output dirctory (default ./output)

-t save p/q value matrix or not
