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
    - [Running the tests](#running-the-tests)
    - [Setting up the database](#setting-up-the-database)
    - [Creates a Database for kmc and ariba](#creates-a-database-for-kmc-and-ariba)
    - [Identify serotype of your input data](#identify-serotype-of-your-input-data)
    - [Summaries the output in one tsv file](#summaries-the-output-in-one-tsv-file)
  - [Output](#output)
  - [Troubleshooting](#troubleshooting)
  - [License](#license)
  - [Citation](#citation)

## Introduction
SeroBA can predict serotypes, by identifying the cps locus, directly from raw whole genome sequencing read data with 98% concordance using a k-mer based method, can process 10,000 samples in just over 1 day using a standard server and can call serotypes at a coverage as low as 10x. SeroBA is implemented in Python3 and is freely available under an open source GPLv3

## Docker Image
Upon each release, a Docker Image is automatically built and pushed to [Docker Hub](https://hub.docker.com/r/sangerbentleygroup/seroba) and [GitHub Packages](https://github.com/sanger-bentley-group/seroba/pkgs/container/seroba)


## Usage
All the following instructions are assuming you are working within a Docker container

### Running the tests
The test can be run from the top level directory:  

```
python3 setup.py test
```

### Setting up the database
SeroBA is packaged with a capsular variant database (CTVdb) which contains references and genetic information for 108 serotypes. It is also possible to add new serotypes by adding the references sequence to the "references.fasta" file in the database folder. Out of the information provided by this database a TSV file is created while using seroba createDBs. You can easily put in additional genetic information for any of these serotypes in the given format.

### Creates a Database for kmc and ariba
```
usage: seroba createDBs  <database dir> <kmer size>

positional arguments:
    database dir     output directory for kmc and ariba Database
    kmer size   kmer_size you want to use for kmc , recommended = 71

Example : 
seroba createDBs my_database/ 71
```
### Identify serotype of your input data
```
usage: seroba runSerotyping [options]  <databases directory> <read1> <read2> <prefix>

    positional arguments:
      database dir         path to database directory
      read1              forward read file
      read2              reverse read file
      prefix             unique prefix

    optional arguments:
      -h, --help         show this help message and exit

    Other options:
      --noclean NOCLEAN  Do not clean up intermediate files (assemblies, ariba
                         report)
      --coverage COVERAGE  threshold for k-mer coverage of the reference sequence (default = 20)                         
```

### Summaries the output in one tsv file
```
usage: seroba summary  <output folder>

positional arguments:
  output folder   directory where the output directories from seroba runSerotyping are stored
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