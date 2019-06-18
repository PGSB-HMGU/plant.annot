# plant.annot
Annotation pipeline for plant genomes by PGSB

# Plant.annot
Repository for the annotation pipeline used @ PGSB. Tested on plant genomes such as Durum, Emmer, T.Urartu, wheat,…

# Requirements

All dependencies should be covered by the provided conda environment 

# Usage

- install anaconda from [anaconda.org](https://anaconda.org)
- install bioconda (see [bioconda](https://bioconda.github.io))
- create plant.annot environment

<code>conda env create --file=plant.annot.yaml</code>

- activate environment

<code>conda activate plant.annot</code>

- install genome threader
- download transposon database PTREP
- download reference proteins from uniprot
- download reference proteins form closely related species
- download additional reference proteins
- edit config.yaml
	- define ISOseq data as described in config.yaml
	- define reference proteins as described in config.yaml
	- define RNAseq data as described in config.yaml
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