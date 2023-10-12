## Nextflow pipeline for processing of GPSeq data

### Requirements
To be able to run this pipeline you need nextflow (version 23.04 or higher) and singularity (tested on version 3.8.6) installed.

The easiest way to install these tools is with conda package manager.

### Running the pipeline
First clone the pipeline (this will create a folder in your current working directory)
`git clone https://github.com/ljwharbers/nextflow-gpseq`
`cd nextflow-gpseq`

You can run the pipeline by typing:
`nextflow run main.nf`
This will run it with default parameters and a test samplesheet and dataset that is included in the repository. If this runs without any issues you can run it in your own dataset.

To do so, you need to create a `samplesheet`. An example samplesheet is located inside the repository and should be a comma separated file that consists of the following columns: `sample,fastq,barcode,condition`.

Following this you can either change the default parameters in the `nextflow.config` file or supply the parameters related to your own dataset in the command you type. I suggest that you change parameters that won't change much between runs in the `nextflow.config` file, while you specify parameters such as `input` and `output` through the command line.

Finally, make sure to also change the `max_memory`, `max_cpus` and `max_time` parameters to make sure you don't go over your systems maximum resources. 

An example command to run this on your own data could be:
`nextflow run main.nf --samplesheet path/to/samplesheet.csv --outdir path/to/results --fasta path/to/reference.fa --bwt2index /path/to/bowtie2/index/folder`

If you do not specify a fasta file and bowtie2 index, you can specify the reference genome you want to use and it will download it from an AWS s3 bucket. For example in the following way:
`nextflow run main.nf --samplesheet path/to/samplesheet.csv --outdir path/to/results --genome GRCh38`

Downloading the fasta file and index might be slow so you can also download the files that you would need through using this tool: https://ewels.github.io/AWS-iGenomes/ Note: you need `aws` tool for this. Once you've downloaded the reference and index files you need you can change the `igenomes_base` parameter in `nextflow.config` and it will take the fasta/index files from there instead of downloading it through nextflow.
