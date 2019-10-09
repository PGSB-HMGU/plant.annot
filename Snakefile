# Created on May 09, 2017
#
# Version 1.0
#
# @author: sven.twardziok@posteo.de


configfile: "ft11.config.yaml"

import csv
from modules import fasta, mygff
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

VERSION = config['cocla']['version']
PREFIX = config['cocla']['prefix']

#####################################################################################################################################
# Hisat2

def hisat2_mapper_input_reads1(wildcards):
    return(config["data"]["rnaseq"][wildcards.dataset][wildcards.sample][1])

def hisat2_mapper_input_reads2(wildcards):
    if len(config["data"]["rnaseq"][wildcards.dataset][wildcards.sample])==1:
        return([])
    if len(config["data"]["rnaseq"][wildcards.dataset][wildcards.sample])==2:
        return(config["data"]["rnaseq"][wildcards.dataset][wildcards.sample][2])

rule hisat2_mapper:
    input:
        reads1 = hisat2_mapper_input_reads1,
        reads2 = hisat2_mapper_input_reads2
    output:
        sam = temp("hisat2/{dataset}/{sample}/{sample}.sam"),
        log = "hisat2/{dataset}/{sample}/{sample}.log"
    params:
        executable = config["executables"]["hisat2"],
        database = config["data"]["hisat2db"],
        arguments = config["hisat2"]["arguments"],
        memory = config["hisat2"]["memory"],
        nodes = config["hisat2"]["nodes"],
        job_name = config['hisat2']['job_name'],
        log = config['hisat2']['log'],
    threads: config["hisat2"]["nodes"]
    resources:
        load = 20 # will run 5 jobs instead of 100
    run:
        if len(config["data"]["rnaseq"][wildcards.dataset][wildcards.sample])==1:
            inreads = ",".join(input.reads1)
            shell("{params.executable} {params.arguments} -p {params.nodes} -x {params.database} -U {inreads} -S {output.sam} 2> {output.log}")
        if len(config["data"]["rnaseq"][wildcards.dataset][wildcards.sample])==2:
            inreads1 = ",".join(input.reads1)
            inreads2 = ",".join(input.reads2)
            shell("{params.executable} {params.arguments} -p {params.nodes} -x {params.database} -1 {inreads1} -2 {inreads2} -S {output.sam} 2> {output.log}")

rule hisat2_sort:
    input:
        sam = "hisat2/{dataset}/{sample}/{sample}.sam"
    output:
        bam = temp("hisat2/{dataset}/{sample}/{sample}_sorted.bam")
    params:
        executable = config["executables"]["samtools"],
        nodes = 1,
        memory= "4G",
        job_name = config['hisat2']['job_name'],
        log = config['hisat2']['log'],
    threads: 8
    resources:
        load = 20
    run:
        shell("{params.executable} sort -@ 8 -o {output.bam} {input.sam}")

rule hisat2_combine:
    input:
        bamfiles=lambda wildcards: ["hisat2/"+ wildcards.dataset + "/" + sample + "/" + sample + "_sorted.bam"
            for sample in config["data"]["rnaseq"][wildcards.dataset].keys()]
    output:
        bam = "hisat2/allsamples_{dataset}.bam"
    params:
        executable = config["executables"]["bamtools"],
        memory="4G",
        nodes = 1,
        job_name = config['hisat2']['job_name'],
        log = config['hisat2']['log'],
    resources:
        load = 1 # will run 100 jobs
    threads: 1
    run:
        shell("{params.executable} merge " + " ".join(["-in " + bamfile for bamfile in input.bamfiles]) + " -out {output.bam}")

#####################################################################################################################################
# Stringtie

rule stringtie:
    input:
        bam = "hisat2/allsamples_{dataset}.bam"
    output:
        gtf = "stringtie/transcripts_{dataset}.gtf"
    params:
        executable = config["executables"]["stringtie"],
        arguments = config["stringtie"]["arguments"],
        memory = config["stringtie"]["memory"],
        nodes = config["stringtie"]["nodes"],
        job_name = config['stringtie']['job_name'],
        log = config['stringtie']['log']
    resources:
        load = 20
    threads: config["stringtie"]["threads"]
    run:
        shell("{params.executable} {params.arguments} -p {params.nodes} -o {output.gtf} {input.bam}")

#####################################################################################################################################
# GMAP

rule gmap_run:
    input:
        fasta=lambda wildcards: config["data"]["longnucl"][wildcards.dataset]
    output:
        gff3 = "gmap/{dataset}_out.gff3"
    params:
        executable = config["executables"]["gmap"],
        dbdir = config["data"]["gmap"]["dbdir"],
        dbname = config["data"]["gmap"]["dbname"],
        arguments = config["gmap"]["arguments"],
        memory = config["gmap"]["memory"],
        nodes = config["gmap"]["nodes"],
        job_name = config['gmap']['job_name'],
        log = config['gmap']['log']
    resources:
        load = 20
    threads: config['gmap']['threads']
    run:
        shell("zcat {input.fasta} | {params.executable} {params.arguments} -f gff3_gene -d {params.dbname} -D {params.dbdir} -t {params.nodes} > {output.gff3}")


rule gmap_combine:
    input:
        gff3files=lambda wildcards: ["gmap/"+ dataset + "_out.gff3"
            for dataset in config["data"]["longnucl"].keys()]
    output:
        gff3file = "gmap/gmap_combined.gff3"
    params:
        nodes = 1,
        memory = "4G",
        job_name = config['gmap']['job_name'],
        log = config['gmap']['log']
    resources:
        load = 1
    threads: 1
    run:
        annos = []
        for gff3file in input.gff3files:
            annos.append(mygff.GeneAnnotation().readGff3PlantAnnot(path=gff3file))
        annos_combined = mygff.GeneAnnotation().combine(annos)
        annos_combined.writeGff3Genes(output.gff3file)

#####################################################################################################################################
# gth

rule gth_splitfasta:
    input:
        fasta = lambda wildcards: config["data"]["refprot"][wildcards.reference]
    output:
        fastas = temp(["gth/{reference}/part_" + str(nbatch) + "/part_" + str(nbatch) + ".fasta"
            for nbatch in range(1, config["gth"]["nbatches"]+1)])
    params:
        nodes = 1,
        memory = "4G",
        job_name = config['gth']['job_name'],
        log = config['gth']['log']
    resources:
        load = 1
    threads: 1
    run:
        splitfasta = fasta.SplitSeqs(sequences=input.fasta, outdir="gth/" + wildcards.reference, nfiles=config["gth"]["nbatches"])

rule gth_refdb:
    input:
        fasta = "gth/{reference}/part_{nbatch}/part_{nbatch}.fasta"
    output:
        dummygenome = temp("gth/{reference}/part_{nbatch}/dummydb.fa"),
        dummygenome_dna_al1 = temp("gth/{reference}/part_{nbatch}/dummydb.fa.dna.al1"),
        dummygenome_dna_des = temp("gth/{reference}/part_{nbatch}/dummydb.fa.dna.des"),
        dummygenome_dna_ois = temp("gth/{reference}/part_{nbatch}/dummydb.fa.dna.ois"),
        dummygenome_dna_prj = temp("gth/{reference}/part_{nbatch}/dummydb.fa.dna.prj"),
        dummygenome_dna_sds = temp("gth/{reference}/part_{nbatch}/dummydb.fa.dna.sds"),
        dummygenome_dna_tis = temp("gth/{reference}/part_{nbatch}/dummydb.fa.dna.tis"),
        dummygenome_md5 = temp("gth/{reference}/part_{nbatch}/dummydb.fa.md5"),
        ref_protein_al1 = temp("gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.al1"),
        ref_protein_bck = temp("gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.bck"),
        ref_protein_bwt = temp("gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.bwt"),
        ref_protein_des = temp("gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.des"),
        ref_protein_lcp = temp("gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.lcp"),
        ref_protein_llv = temp("gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.llv"),
        ref_protein_ois = temp("gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.ois"),
        ref_protein_prj = temp("gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.prj"),
        ref_protein_sds = temp("gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.sds"),
        ref_protein_ssp = temp("gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.ssp"),
        ref_protein_sti1 = temp("gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.sti1"),
        ref_protein_suf = temp("gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.suf"),
        ref_protein_tis = temp("gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.tis"),
        ref_md5 = temp("gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.md5")
    params:
        executable=config["executables"]["gth"],
        nodes = 1,
        memory = "4G",
        job_name = config['gth']['job_name'],
        log = config['gth']['log']
    resources:
        load = 1
    threads: 1
    run:
        dummydb =  output.dummygenome
        with open(dummydb, "w") as dummyfile:
            dummyfile.write(">dummy_db\n")
            dummyfile.write("ACGTCACGACGAGAATAGTGTAAACAAAACTGCTGTCGGCGGAAGCGCCAAAGGAGTCTGTGAATTCTTATTCCCGAATAACATCCGTCTCCGTGCGGGAAAATCACCGACGCCGTTTTATAGAAG\n")
        shell(params.executable + " -createindicesonly -genomic " + dummydb + " -protein " + input.fasta + " -species rice")


rule gth_run:
    input:
        ref="gth/{reference}/part_{nbatch}/part_{nbatch}.fasta",
        ref_protein_al1 = "gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.al1",
        ref_protein_bck = "gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.bck",
        ref_protein_bwt = "gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.bwt",
        ref_protein_des = "gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.des",
        ref_protein_lcp = "gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.lcp",
        ref_protein_llv = "gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.llv",
        ref_protein_ois = "gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.ois",
        ref_protein_prj = "gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.prj",
        ref_protein_sds = "gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.sds",
        ref_protein_ssp = "gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.ssp",
        ref_protein_sti1 = "gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.sti1",
        ref_protein_suf = "gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.suf",
        ref_protein_tis = "gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.protein.tis",
        ref_md5 = "gth/{reference}/part_{nbatch}/part_{nbatch}.fasta.md5",
        database=lambda wildcards: config["data"]["gth"][wildcards.database]
    output:
        gff3 = "gth/{reference}/part_{nbatch}/{database}/gth_genes.gff3"
    params:
        executable = config["executables"]["gth"],
        arguments = config["gth"]["arguments"],
        memory = config["gth"]["memory"],
        nodes = config["gth"]["nodes"],
        job_name = "gthing",
        log = config['gth']['log']
    resources:
        load = 1
    threads: 1
    run:
        shell(params.executable + " -skipindexcheck -genomic " + input.database + " -protein " + input.ref +
            " -skipalignmentout {params.arguments} -gff3out -force -o " + output.gff3)

rule gth_combine_dbs:
    input:
        gff3files = lambda wildcards: ["gth/"+ wildcards.reference + "/part_" +  wildcards.nbatch + "/" + database + "/gth_genes.gff3"
            for database in config["data"]["gth"].keys()]
    output:
        gff3file = "gth/{reference}/part_{nbatch}/gth_genes_combined.gff3"
    params:
        memory = "4G",
        nodes = 1,
        job_name = config['gth']['job_name'],
        log = config['gth']['log']
    resources:
        load = 1
    threads: 1
    run:
        annos = []
        for gff3file in input.gff3files:
            annos.append(mygff.GeneAnnotation().readGff3PlantAnnot(path=gff3file))
        annos_combined = mygff.GeneAnnotation().combine(annos)
        annos_combined.writeGff3Genes(output.gff3file)

rule gth_combine_dataset:
    input:
        gff3files = lambda wildcards: ["gth/"+ wildcards.reference + "/part_" +  str(nbatch) + "/gth_genes_combined.gff3"
            for nbatch in range(1, config["gth"]["nbatches"]+1)]
    output:
        gff3file = "gth/genes_{reference}.gff3"
    params:
        memory = "4G",
        nodes = 1,
        job_name = config['gth']['job_name'],
        log = config['gth']['log']
    resources:
        load = 1
    threads: 1
    run:
        annos = []
        for gff3file in input.gff3files:
            annos.append(mygff.GeneAnnotation().readGff3PlantAnnot(path=gff3file))
        annos_combined = mygff.GeneAnnotation().combine(annos)
        annos_combined.writeGff3Genes(output.gff3file)

#####################################################################################################################################
# Transdecoder

def transdecoder_cuffcompare_input(wildcards):
    tmpgffs = []
    if "refprot" in config["data"].keys():
        for reference in config["data"]["refprot"]:
            tmpgffs.append("gth/genes_" + reference + ".gff3")
    if "rnaseq" in config["data"].keys():
        for dataset in config["data"]["rnaseq"]:
            tmpgffs.append("stringtie/transcripts_" + dataset + ".gtf")
    if "longnucl" in  config["data"].keys():
        tmpgffs.append("gmap/gmap_combined.gff3")
    return(tmpgffs)

rule transdecoder_cuffcompare:
    input:
        gffs = transdecoder_cuffcompare_input
    output:
        gtf=temp("evidences.combined.gtf"),
        loci=temp("evidences.loci"),
        stats=temp("evidences.stats"),
        tracking=temp("evidences.tracking")
    params:
        executable=config["executables"]["cuffcompare"],
        nodes = 1,
        memory = "4G",
        job_name = "cuffcompare",
        log = config['transdecoder']['log']
    resources:
        load = 1
    threads: 1
    run:
        tmpgffs = " ".join(input.gffs)
        shell("{params.executable} -o evidences " + tmpgffs)

rule transdecoder_stringtie:
    input:
        gtf = "evidences.combined.gtf"
    output:
        gtf = temp("evidences.combined.stringtie_merged.gtf")
    params:
        executable=config["executables"]["stringtie"],
        nodes = config["transdecoder"]["stringtie"]["nodes"],
        memory = config["transdecoder"]["stringtie"]["memory"],
        job_name = "stringtie",
        log = config['transdecoder']['log']
    resources:
        load = 1
    threads: 1
    run:
        shell("{params.executable} --merge -o {output.gtf} {input.gtf}")

rule transdecoder_extract:
    input:
        gtf="evidences.combined.stringtie_merged.gtf"
    output:
        fasta= "transcripts.fasta"
    params:
        executable=config["executables"]["transdecoder"]["extract"],
        genome=config["data"]["genome"],
        nodes = 1,
        memory = "4G",
        job_name = config['transdecoder']['job_name'],
        log = config['transdecoder']['log']
    resources:
        load = 1
    threads: 1
    run:
        shell("{params.executable} {input.gtf} {params.genome} > {output.fasta}")

rule transdecoder_longorfs:
    input:
        fasta="transcripts.fasta"
    output:
        pep="transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    params:
        executable=config["executables"]["transdecoder"]["longorfs"],
        genome=config["data"]["genome"],
        nodes = 1,
        memory = "4G",
        job_name = "longorfs",
        log = config['transdecoder']['log']
    resources:
        load = 1
    threads: 1
    run:
        shell("{params.executable} -t {input.fasta} -p 0")

rule transdecoder_splitfasta:
    input:
        fasta="transcripts.fasta.transdecoder_dir/longest_orfs.pep"
    output:
        fastas=temp(["transcripts.fasta.transdecoder_dir/batches/part_" + str(nbatch) + "/part_" + str(nbatch) + ".fasta"
            for nbatch in range(1, config["transdecoder"]["nbatches"]+1)])
    params:
        nodes = 1,
        memory = "4G",
        job_name = config['transdecoder']['job_name'],
        log = config['transdecoder']['log']
    resources:
        load = 1
    threads: 1
    run:
        splitfasta = fasta.SplitSeqs(sequences=input.fasta, outdir="transcripts.fasta.transdecoder_dir/batches" , nfiles=config["transdecoder"]["nbatches"])

rule transdecoder_blast:
    input:
        fasta="transcripts.fasta.transdecoder_dir/batches/part_{nbatch}/part_{nbatch}.fasta"
    output:
        blp=temp("transcripts.fasta.transdecoder_dir/batches/part_{nbatch}/part_{nbatch}.blp")
    params:
        executable = config["executables"]["blastp"],
        database = config["data"]["transdecoder"]["blastp"],
        nodes = config["transdecoder"]["blastp"]["nodes"],
        memory = config["transdecoder"]["blastp"]["memory"],
        job_name = "blasting",
        log = config['transdecoder']['log']
    resources:
        load = 1
    threads: 1
    run:
        shell(params.executable + " -max_target_seqs 1 -evalue 1e-05 -db {params.database} -query {input.fasta} -out {output.blp} -outfmt 6")

rule transdecoder_blast_combine:
    input:
        blps=lambda wildcards: ["transcripts.fasta.transdecoder_dir/batches/part_" + str(nbatch) + "/part_" + str(nbatch) + ".blp"
            for nbatch in range(1, config["transdecoder"]["nbatches"]+1)]
    output:
        blp="transcripts.fasta.transdecoder_dir/longest_orfs.pep_blastresults.blp"
    params:
        nodes = 1,
        memory = "4G",
        job_name = config['transdecoder']['job_name'],
        log = config['transdecoder']['log']
    resources:
        load = 1
    threads: 1
    run:
        shell("touch {output.blp}")
        for blp in input.blps:
             shell("cat " + blp + " >> {output.blp}")

rule transdecoder_hmmscan:
    input:
        fasta="transcripts.fasta.transdecoder_dir/batches/part_{nbatch}/part_{nbatch}.fasta"
    output:
        domtblout=temp("transcripts.fasta.transdecoder_dir/batches/part_{nbatch}/part_{nbatch}.domtblout")
    params:
        executable = config["executables"]["hmmscan"],
        pfamhmm = config["data"]["transdecoder"]["pfamhmm"],
        nodes = config["transdecoder"]["hmmscan"]["nodes"],
        memory = config["transdecoder"]["hmmscan"]["memory"],
        job_name = "hmmscanning",
        log = config['transdecoder']['log']
    resources:
        MB = 2000,
        load = 1
    threads: 2
    run:
        shell(params.executable + "  --domtblout {output.domtblout} {params.pfamhmm} {input.fasta}")

rule transdecoder_hmmscan_combine:
    input:
        domtblout=lambda wildcards: ["transcripts.fasta.transdecoder_dir/batches/part_" + str(nbatch) + "/part_" + str(nbatch) + ".domtblout"
            for nbatch in range(1, config["transdecoder"]["nbatches"]+1)]
    output:
        domtblout="transcripts.fasta.transdecoder_dir/longest_orfs.pep_hmmscan.domtblout"
    params:
        nodes = 1,
        memory = "4G",
        job_name = config['transdecoder']['job_name'],
        log = config['transdecoder']['log']
    resources:
        load = 1,
        MB = 2000
    threads: 1
    run:
        shell("touch {output.domtblout}")
        for domtblout in input.domtblout:
             shell("grep -v \"#\" " + domtblout + " >> {output.domtblout}")

rule transdecoder_predict:
    input:
        fasta = "transcripts.fasta",
        blp = "transcripts.fasta.transdecoder_dir/longest_orfs.pep_blastresults.blp",
        domtblout = "transcripts.fasta.transdecoder_dir/longest_orfs.pep_hmmscan.domtblout"
    output:
        gff3 = "transcripts.fasta.transdecoder.gff3"
    params:
        executable=config["executables"]["transdecoder"]["predict"],
        nodes = config["transdecoder"]["predict"]["nodes"],
        memory = config["transdecoder"]["predict"]["memory"],
        job_name = "predicting",
        log = config['transdecoder']['log']
    resources:
        load = 1
    threads: 1
    run:
        shell("{params.executable} -t {input.fasta} --retain_pfam_hits {input.domtblout} --retain_blastp_hits {input.blp} --cpu {params.nodes}")

rule transdecoder_convert:
    input:
        fasta = "transcripts.fasta",
        gff3 = "transcripts.fasta.transdecoder.gff3",
        gtf="evidences.combined.stringtie_merged.gtf"
    output:
        gff3tmp = "transcripts.gff3",
        gff3 = "transcripts.genes.gff3"
    params:
        executable_gff3=config["executables"]["transdecoder"]["convertgff3"],
        executable_genome=config["executables"]["transdecoder"]["convertgenome"],
        nodes = config["transdecoder"]["convert"]["nodes"],
        memory = config["transdecoder"]["convert"]["memory"],
        job_name = "converting",
        log = config['transdecoder']['log']
    resources:
        load = 1
    threads: 1
    run:
        shell("{params.executable_gff3} {input.gtf} > {output.gff3tmp}")
        shell("{params.executable_genome} {input.gff3} {output.gff3tmp} {input.fasta} > {output.gff3}")

#####################################################################################################################################
# combine transdecoder results with genomethreader results and remove equal protein sequences

def candidates_combine_input(wildcards):
    tmpgffs = ["transcripts.genes.gff3"]
    if "refprot" in config["data"].keys():
        for reference in config["data"]["refprot"]:
            tmpgffs.append("gth/genes_" + reference + ".gff3")
    return(tmpgffs)

rule candidates_combine:
    input:
        gff3files = candidates_combine_input
    output:
        gff3 = "transcripts.genes.combined.gff3",
        cds = "transcripts.genes.combined_cds.fasta",
        protein = "transcripts.genes.combined_protein.fasta"
    params:
        executable = config["executables"]["gffread"],
        genome = config["data"]["genome"],
        memory = "8G",
        nodes = 1,
        job_name = "combining",
        log = "combine.log"
    resources:
        load = 1
    threads: 1
    run:
        print("read gff3 files and create annotations")
        annos = []
        for gff3file in input.gff3files:
            print("    %s" % (gff3file))
            annos.append(mygff.GeneAnnotation().readGff3PlantAnnot(path=gff3file))
        print("Combine annotations")
        annos_combined = mygff.GeneAnnotation().combine(annos)
        print("Calculate new geneids")
        anno_newgeneids = annos_combined.recalcGeneids()
        print("Write annotation")
        anno_newgeneids.writeGff3Genes(output.gff3)
        print("Write coding sequences")
        shell("{params.executable} {output.gff3} -g {params.genome} -x {output.cds}")
        print("Write proteins")
        with open(output.cds, "r") as infile:
            with open(output.protein, "w") as outfile:
                for record in SeqIO.parse(infile, "fasta"):
                    if len(str(record.seq)) % 3 == 0:
                        tmprecord = SeqRecord(record.seq.translate(), id=record.id, description="")
                        SeqIO.write(tmprecord, outfile, "fasta")

rule candidates_collapse:
    input:
        gff3 = "transcripts.genes.combined.gff3",
    output:
        gff3 = "transcripts.genes.combined.collapsed.gff3",
        cds = "transcripts.genes.combined.collapsed_cds.fasta",
        protein = "transcripts.genes.combined.collapsed_protein.fasta"
    params:
        executable = config["executables"]["gffread"],
        genome = config["data"]["genome"],
        nodes = 1,
        memory = "16G",
        job_name = "collapsing",
        log = "collapsing.log"
    resources:
        load = 1
    threads: 1
    run:
        print("Read annotation")
        anno = mygff.GeneAnnotation().readGff3PlantAnnot(path=input.gff3)
        print("Collapse models")
        anno_collapsed = anno.collapseMrnas()
        print("Write annotation")
        anno_collapsed.writeGff3Genes(output.gff3)
        print("Write coding sequences")
        shell("{params.executable} {output.gff3} -g {params.genome} -x {output.cds}")
        print("Write proteins")
        with open(output.cds, "r") as infile:
            with open(output.protein, "w") as outfile:
                for record in SeqIO.parse(infile, "fasta"):
                    if len(str(record.seq)) % 3 == 0:
                        tmprecord = SeqRecord(record.seq.translate(), id=record.id, description="")
                        SeqIO.write(tmprecord, outfile, "fasta")

rule candidates_stats:
    input:
        gff3 = "transcripts.genes.combined.collapsed.gff3",
        cds = "transcripts.genes.combined.collapsed_cds.fasta",
    output:
        genes = "transcripts.genes.combined.collapsed_genes_stats.tab",
        transcripts = "transcripts.genes.combined.collapsed_transcripts_stats.tab",
        cds = "transcripts.genes.combined.collapsed_cds_stats.tab"
    params:
        nodes = 1,
        memory = "16G",
        job_name = "stats",
        log = "stats.log"
    resources:
        load = 1
    threads: 1
    run:
        print("Write CDS stats")
        fasta.PrintCdsStats(infasta=input.cds, outstats=output.cds)
        print("Read annotation")
        anno = mygff.GeneAnnotation().readGff3PlantAnnot(path=input.gff3)
        print("Write genes stats")
        anno.printGeneStats(path=output.genes)
        print("Write transcripts stats")
        anno.printTranscriptsStats(path=output.transcripts)

#####################################################################################################################################
# confidence classification

rule cocla_splitfasta:
    input:
        fasta="transcripts.genes.combined.collapsed_protein.fasta"
    output:
        fastas=["cocla/{database}/batches/part_" + str(nbatch) + "/part_" + str(nbatch) + ".fasta"
            for nbatch in range(1, config["cocla"]["nbatches"]+1)]
    params:
        nodes = 1,
        memory = "20G",
        job_name = "chunking",
        log = "cocla.log"
    resources:
        load = 1
    threads: 1
    run:
        splitfasta = fasta.SplitSeqs(sequences=input.fasta, outdir="cocla/"+ wildcards.database +"/batches" , nfiles=config["cocla"]["nbatches"])

rule cocla_blast:
    input:
        fasta = "cocla/{database}/batches/part_{nbatch}/part_{nbatch}.fasta",
        database = lambda wildcards: config["data"]["cocla"][wildcards.database]
    output:
        blp="cocla/{database}/batches/part_{nbatch}/part_{nbatch}.blp"
    params:
        executable = config["executables"]["blastp"],
        nodes = config["cocla"]["nodes"],
        memory = config["cocla"]["memory"],
        evalue = config["cocla"]["evalue"],
        job_name = "cocla_blast",
        log = "cocla.log"
    resources:
        load = 1
    threads: 1
    run:
        shell(params.executable + " -evalue {params.evalue}  -num_threads {params.nodes} -subject {input.database} -query {input.fasta} -out {output.blp} -outfmt \"6 qseqid sseqid pident qlen qstart qend slen sstart send evalue bitscore\"")

rule cocla_blast_combine:
    input:
        blps=lambda wildcards: ["cocla/" + wildcards.database + "/batches/part_" + str(nbatch) + "/part_" + str(nbatch) + ".blp"
            for nbatch in range(1, config["cocla"]["nbatches"]+1)]
    output:
        blp="cocla/results_{database}.blp"
    params:
        nodes = 1,
        memory = "16G",
        job_name = "combining",
        log = "cocla.log"
    resources:
        load = 1
    threads: 1
    run:
        shell("touch {output.blp}")
        for blp in input.blps:
             shell("cat " + blp + " >> {output.blp}")


rule cocla_all:
    input:
        transcripts = PREFIX+"_"+VERSION+"transcripts.genes.combined.collapsed_transcripts_stats.tab",
        cds = PREFIX+"_"+VERSION+"transcripts.genes.combined.collapsed_cds_stats.tab",
        trep = "cocla/results_trep.blp",
        unimag = "cocla/results_unimag.blp",
        unipoa = "cocla/results_unipoa.blp"

rule cocla_report:
    input:
        transcripts = "transcripts.genes.combined.collapsed_transcripts_stats.tab",
        cds = "transcripts.genes.combined.collapsed_cds_stats.tab",
        trep = "cocla/results_trep.blp",
        unimag = "cocla/results_unimag.blp",
        unipoa = "cocla/results_unipoa.blp"
    output:
        report = "confidence_classification_"+VERSION+".html"
    params:
        nodes = 1,
        memory = "16G",
        job_name = "reporting",
        log = "cocla.log",
        unipoa_threshold = config['cocla']['unipoa_threshold'],
        unimag_threshold = config['cocla']['unimag_threshold'],
        repeat_threshold = config['cocla']['repeat_threshold'],
    resources:
        load = 1
    threads: 1
    script:
        "confidence_classification.Rmd"

#####################################################################################################################################
# get final files

rule final_files:
    input:
        gff3 = "transcripts.genes.combined.collapsed.gff3",
        report = "confidence_classification_"+VERSION+".html"
    output:
        gff3 = PREFIX+"_"+VERSION+".transcripts.genes.combined.collapsed_cocla.gff3",
        gff3HC = PREFIX+"_"+VERSION+".transcripts.genes.combined.collapsed_cocla_HC.gff3",
        gff3LC = PREFIX+"_"+VERSION+".transcripts.genes.combined.collapsed_cocla_LC.gff3",
        cds = PREFIX+"_"+VERSION+".transcripts.genes.combined.collapsed_cocla_cds.fasta",
        cdsHC = PREFIX+"_"+VERSION+".transcripts.genes.combined.collapsed_cocla_HC_cds.fasta",
        cdsLC = PREFIX+"_"+VERSION+".transcripts.genes.combined.collapsed_cocla_LC_cds.fasta",
        protein = PREFIX+"_"+VERSION+".transcripts.genes.combined.collapsed_cocla_cds_proteins.fasta",
        proteinHC = PREFIX+"_"+VERSION+".transcripts.genes.combined.collapsed_cocla_HC_cds_proteins.fasta",
        proteinLC = PREFIX+"_"+VERSION+".transcripts.genes.combined.collapsed_cocla_LC_cds_proteins.fasta"
    params:
        executable = config["executables"]["gffread"],
        genome = config["data"]["genome"],
        nodes = 1,
        memory = "16G",
        job_name = "final_files",
        log = "cocla.log"
    resources:
        load = 1
    threads: 1
    run:
        conflcass = {}
        with open("confidence_classification_"+VERSION+".tab", "r") as csvfile:
            reader = csv.DictReader(csvfile, delimiter=",")
            for line in reader:
                conflcass[line["id"]] = {}
                conflcass[line["id"]]["primconf"] = line["primconf"]
                conflcass[line["id"]]["secconf"] = line["secconf"]
        anno = mygff.GeneAnnotation().readGff3PlantAnnot(path=input.gff3)
        newfeatures = []
        newgenes = {}
        newmrnas = {}
        newseqids = {}
        for geneid in anno.genes:
            gene = anno.genes[geneid]
            tmpclass = "none"
            tmphcmrnas = []
            tmplcmrnas = []
            for mrna in gene.features:
                mrna.primary_confidence_class = conflcass[mrna.identifier]["primconf"]
                mrna.secondary_condidence_class = conflcass[mrna.identifier]["secconf"]
                if conflcass[mrna.identifier]["primconf"] == "HC":
                    tmpclass = "HC"
                    tmphcmrnas.append(mrna)
                if tmpclass != "HC" and conflcass[mrna.identifier]["primconf"] == "LC":
                    tmpclass = "LC"
                    tmplcmrnas.append(mrna)
            if tmpclass=="HC":
                gene.features = tmphcmrnas
                gene.primary_confidence_class = "HC"
                newfeatures.append(gene)
                newgenes[gene.identifier] = gene
                for mrna in tmphcmrnas:
                    newfeatures.append(mrna)
                    newmrnas[mrna.identifier] = mrna
                    newfeatures += mrna.features
            if tmpclass=="LC":
                gene.mrnas = tmplcmrnas
                gene.primary_confidence_class = "LC"
                newfeatures.append(gene)
                newgenes[gene.identifier] = gene
                for mrna in tmplcmrnas:
                    newfeatures.append(mrna)
                    newmrnas[mrna.identifier] = mrna
                    newfeatures += mrna.features
        anno.features = sorted(newfeatures)
        anno.genes = newgenes
        anno.mrnas = newmrnas
        print("Write new gff3 files")
        anno_newgenes = anno.recalcGeneids()
        anno_newgenes.writeGff3Genes(output.gff3)
        anno_hc = anno.getHcGff3Genes()
        anno_hc.writeGff3Genes(output.gff3HC)
        anno_lc = anno.getLcGff3Genes()
        anno_lc.writeGff3Genes(output.gff3LC)
        print("Write coding sequences")
        shell("{params.executable} {output.gff3} -g {params.genome} -x {output.cds}")
        shell("{params.executable} {output.gff3HC} -g {params.genome} -x {output.cdsHC}")
        shell("{params.executable} {output.gff3LC} -g {params.genome} -x {output.cdsLC}")
        print("Write proteins sequences")
        with open(output.cds, "r") as infile:
            with open(output.protein, "w") as outfile:
                for record in SeqIO.parse(infile, "fasta"):
                    if len(str(record.seq)) % 3 == 0:
                        tmprecord = SeqRecord(record.seq.translate(), id=record.id, description="")
                        SeqIO.write(tmprecord, outfile, "fasta")
        with open(output.cdsHC, "r") as infile:
            with open(output.proteinHC, "w") as outfile:
                for record in SeqIO.parse(infile, "fasta"):
                    if len(str(record.seq)) % 3 == 0:
                        tmprecord = SeqRecord(record.seq.translate(), id=record.id, description="")
                        SeqIO.write(tmprecord, outfile, "fasta")
        with open(output.cdsLC, "r") as infile:
            with open(output.proteinLC, "w") as outfile:
                for record in SeqIO.parse(infile, "fasta"):
                    if len(str(record.seq)) % 3 == 0:
                        tmprecord = SeqRecord(record.seq.translate(), id=record.id, description="")
                        SeqIO.write(tmprecord, outfile, "fasta")
