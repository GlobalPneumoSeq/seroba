# SeroBA
This is a fork of [Wellcome Sanger Institute Pathogen Informatics' SeroBA](https://github.com/sanger-pathogens/seroba). As the original SeroBA is no longer maintained, this fork mainly aims to integrate bug fixes and database updates to provide a stable, updated, and containerised version of SeroBA. 

## About 
SeroBA is a k-mer based Pipeline to identify the Serotype from Illumina NGS reads for given references.

## Contents
- [SeroBA](#seroba)
  - [About](#about)
  - [Contents](#contents)
  - [Introduction](#introduction)
  - [Docker Image](#docker-image)
  - [Usage](#usage)
    - [Running with docker](#running-with-docker)
    - [Running with singularity](#running-with-singularity)
    - [Running the tests](#running-the-tests)
    - [Running on multiple samples](#run-on-multiple-samples)
    - [Summarise the output in one csv file](#summarise-the-output-in-one-csv-file)
  - [Output](#output)
  - [Troubleshooting](#troubleshooting)
  - [License](#license)
  - [Citation](#citation)

## Introduction
SeroBA can predict serotypes, by identifying the cps locus, directly from raw whole genome sequencing read data with 98% concordance using a k-mer based method, can process 10,000 samples in just over 1 day using a standard server and can call serotypes at a coverage as low as 10x. SeroBA is implemented in Python3 and is freely available under an open source GPLv3

## Docker Image
Upon each release, a Docker Image is automatically built and pushed to [Docker Hub](https://hub.docker.com/r/sangerbentleygroup/seroba) and [GitHub Packages](https://github.com/sanger-bentley-group/seroba/pkgs/container/seroba)


## Usage

### Running with docker
To run serotyping with docker using the pre-built docker image which contains the database, run a command like below. Replace the placeholder values (/path/to/reads, read1_file_name, read2_file_name, output_folder_prefix) in the command.

```
docker run --rm -it -u $(id -u):$(id -g) -v /path/to/reads:/data sangerbentleygroup/seroba seroba runSerotyping /seroba/database /data/read1_file_name /data/read2_file_name /data/output_folder_prefix
```

### Running with singularity
To run serotyping with singularity using the pre-built docker image which contains the database, run a command like below. Replace the placeholder values (/path/to/reads, read1_file_name, read2_file_name, output_folder_prefix) in the command.

```
singularity exec --bind /path/to/reads:/data docker://sangerbentleygroup/seroba seroba runSerotyping /seroba/database /data/read1_file_name /data/read2_file_name /data/output_folder_prefix
```

### Run on multiple samples
1. Place all your files into a single directory and `cd` into that directory
2. Run seroBA on all samples using a for loop

Run with docker:
```
for READ1 in *1.fastq.gz; do SAMPLE=${READ1%_1.fastq.gz}; docker run --rm -it -u $(id -u):$(id -g) -v $PWD:/data sangerbentleygroup/seroba seroba runSerotyping /seroba/database /data/${SAMPLE}_1.fastq.gz /data/${SAMPLE}_2.fastq.gz /data/${SAMPLE}_RESULT; done
```

Run with singularity:
```
for READ1 in *1.fastq.gz; do SAMPLE=${READ1%_1.fastq.gz}; singularity exec --bind $PWD:/data docker://sangerbentleygroup/seroba seroba runSerotyping /seroba/database /data/${SAMPLE}_1.fastq.gz /data/${SAMPLE}_2.fastq.gz /data/${SAMPLE}_RESULT; done
```

### Summarise the output in one CSV file
To summarise all the output into one CSV file, you can use the `seroba summary` command, run from the folder where your seroBA results are stored:

Run with docker:
```
docker run --rm -it -u $(id -u):$(id -g) -v $PWD:/data sangerbentleygroup/seroba seroba summary /data/
```

Run with singularity:
```
singularity exec --bind $PWD:/data docker://sangerbentleygroup/seroba seroba summary /data/
```

The summary file will be available as `summary.csv` in the directory.

### Running the tests
To run the tests using docker, run the below command:
```
docker run --workdir /seroba -it --rm sangerbentleygroup/seroba python3 setup.py test
```


## Output
In the folder 'prefix' you will find a file named `pred.csv` including your predicted serotype and genetic variant as well as a file called detailed_serogroup_info.txt including information about SNP, genes, and alleles that are found in your reads. After the use of `seroba summary` a csv file called `summary.csv` is created that consists of four columns (Sample,Serotype,Genetic_Variant,Contamination_Status). Serotypes that do not match any reference are marked as "untypable".

__detailed_serogroup_info example:__
```
Predicted Serotype:       23F
Serotype predicted by ariba:    23F
assembly from ariba has an identity of:   99.77    with this serotype

Serotype       Genetic Variant
23F            allele  wchA
```
In the detailed information you can see the finally predicted serotype as well as the serotypes that had the closest reference in that specific serogroup according to ARIBA. Furthermore you can see the sequence identity between the sequence assembly and the reference sequence.  

## Troubleshooting
* Case 1:
	* SeroBA predicts 'untypable'. An 'untypable' prediction can either be a
real 'untypable' strain or can be caused by different problems. Possible problems are:
bad quality of your input data, submission of a wrong species or to low coverage
of your sequenced reads. Please check your data again and run a quality control.

* Case 2:
	* 	Low alignment identity in the 'detailed_serogroup_info' file. This can
be a hint for a mosaic serotpye.
	* Possible solution: perform a blast search on the whole genome assembly

* Case 3:
	* The fourth column in the summary.csv indicates "contamination". This means that
    at least one heterozygous SNP was detected in the read data with at least
    10% of the mapped reads at the specific position supporting the SNP.
	* Possible solution: please check the quality of your data and have a look
    for contamination within your reads

## License
SeroBA is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/seroba/blob/master/LICENSE)

## Citation
__SeroBA: rapid high-throughput serotyping of Streptococcus pneumoniae from whole genome sequence data__  
Epping L, van Tonder, AJ, Gladstone RA, GPS Consortium, Bentley SD, Page AJ, Keane JA, Microbial Genomics 2018, doi: [10.1099/mgen.0.000186](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000186)