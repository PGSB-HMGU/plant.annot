# Plant.annot
Repository for the annotation pipeline used @ PGSB. Tested on plant genomes.


# Requirements

Most dependencies should be covered by the provided conda environment 

# Usage

- install anaconda from [anaconda.org](https://anaconda.org)
- install bioconda (see [bioconda](https://bioconda.github.io))
- clone repository

<code>git clone https://github.com/PGSB-HMGU/plant.annot.git</code>

- cd into plant.annot directory

<code>cd plant.annot</code>



- create plant.annot environment using the provided yaml file

<code>conda env create --file=plant.annot.yaml</code>

- activate environment

<code>conda activate plant.annot</code>

- download and install Transdecoder 3.0.1 from [Transdecoder](https://github.com/TransDecoder/TransDecoder/wiki)
- download transposon database PTREP from [Hypothetical TREP protein sequences](https://botserv2.uzh.ch/kelldata/trep-db/downloadFiles.html)
- download reference proteins from uniprot
- download reference proteins form closely related species
- edit config.yaml
    - define ISOseq data as described in config.yaml
    - define reference proteins as described in config.yaml
    - define RNAseq data as described in config.yaml
    - review executable section in config.yaml 
    - build hisat2 index
    - build gmap index
    - split large chromosomes into single files to speed up gth
    - 
- perform a dry run

<code>snakemake final_files -np</code>

- run locally

<code>snakemake final_files --cores *cores*</code>

- run on a cluster

<code>snakemake fina_files --cluster "qsub"</code>