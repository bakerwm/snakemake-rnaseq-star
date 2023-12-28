# Prepare reference for RNAseq pipeline

## Output files/directories in this snakemake
# fasta:         resources/{ref_build}/genome.fasta
# gtf:           resources/{ref_build}/genome.gtf
# bowtie2_index: resources/{ref_build}/genome.fasta
# bwa_index:     resources/{ref_build}/genome.fasta
# star_index:    resources/{ref_build}/star_genome
# hisat2_index:  resources/{ref_build}/genome
# salmon_index:  resources/{ref_build}/salmon_genome
# decoys:        resources/{ref_build}/genome.decoys.txt
# cdna:          resources/{ref_build}/genome.cdna.fasta
# gentrome:      resources/{ref_build}/genome.transcriptome.fasta
# exons:         resources/{ref_build}/genome.exons.bed
# splice_sites:  resources/{ref_build}/genome.splice_sites.bed

###############################################################################
## local variables
REF_RELEASE = config["ref"]["release"]
REF_BUILD = config["ref"]["build"]  # genome build, for sub_dir of genome
REF_SPECIES = config["ref"]["species"]  # species name
THREADS = 8
## Wrappers from Github
# WRAPPER_PATH = "v3.3.0" # wrappers from Github
## Wrappers from local machine, public directory, version 3.3.2, 2023-12-26
WRAPPER_PATH = "file:///data/biosoft/snakemake/snakemake-wrappers"


###############################################################################
# Choose ref data from config or build in this pipeline
# required:
# - fasta
# - gtf
# - star_index (!!! caution: STAR version conflic)
# - salmon_index
def get_ref_fasta():
    cfg_fa = config["index"]["fasta"]
    out = f"resources/{REF_BUILD}/genome.fasta"
    if isinstance(cfg_fa, str):
        if Path(cfg_fa).exists():
            out = cfg_fa
    return out


def get_ref_gtf():
    cfg_gtf = config["index"]["gtf"]
    out = f"resources/{REF_BUILD}/genome.gtf"
    if isinstance(cfg_gtf, str):
        if Path(cfg_gtf).exists():
            out = cfg_gtf
    return out


def get_ref_star_index():
    cfg_star_index = config["index"]["star_index"]
    out = f"resources/{REF_BUILD}/star_genome"
    # STAR --runMode genomeLoad --genomeDir /path/to/your/genome/index
    if isinstance(cfg_star_index, str):
        out = cfg_star_index
    return out


def get_ref_salmon_index():
    cfg_salmon_index = config["index"]["salmon_index"]
    out = f"resources/{REF_BUILD}/salmon_genome"
    if isinstance(cfg_salmon_index, str):
        out = cfg_salmon_index
    return out


###############################################################################
# Download and build reference fasta and index
rule get_genome:
    output:
        f"resources/{REF_BUILD}/genome.fasta",
    log:
        f"logs/genome/get-genome.{REF_BUILD}.log",
    params:
        species=REF_SPECIES,
        datatype="dna",
        build=REF_BUILD,
        release=REF_RELEASE,
    cache: True
    wrapper:
        f"{WRAPPER_PATH}/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        f"resources/{REF_BUILD}/genome.gtf",
    params:
        species=REF_SPECIES,
        fmt="gtf",
        build=REF_BUILD,
        release=REF_RELEASE,
        flavor="",
    cache: True
    log:
        f"logs/genome/get_annotation.{REF_BUILD}.log",
    wrapper:
        f"{WRAPPER_PATH}/bio/reference/ensembl-annotation"


rule get_transcriptome:
    output:
        f"resources/{REF_BUILD}/genome.cdna.fasta",
    params:
        species=REF_SPECIES,
        build=REF_BUILD,
        release=REF_RELEASE,
        datatype="cdna",
    cache: True
    log:
        f"logs/genome/get_transcriptome.{REF_BUILD}.log",
    wrapper:
        f"{WRAPPER_PATH}/bio/reference/ensembl-sequence"


rule genome_faidx:
    input:
        f"resources/{REF_BUILD}/genome.fasta",
    output:
        f"resources/{REF_BUILD}/genome.fasta.fai",
    log:
        f"logs/genome/get-faidx.{REF_BUILD}.log",
    cache: True
    wrapper:
        f"{WRAPPER_PATH}/bio/samtools/faidx"


rule bwa_index:
    input:
        f"resources/{REF_BUILD}/genome.fasta",
    output:
        multiext(
            f"resources/{REF_BUILD}/genome",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        f"logs/genome/bwa_index.{REF_BUILD}.log",
    threads: 8
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        f"{WRAPPER_PATH}/bio/bwa/index"


rule star_index:
    input:
        fasta=f"resources/{REF_BUILD}/genome.fasta",
        annotation=f"resources/{REF_BUILD}/genome.gtf",
    output:
        directory(f"resources/{REF_BUILD}/star_genome"),
    threads: 24
    params:
        # extra=f"--sjdbGTFfile resources/{REF_BUILD}/genome.gtf --sjdbOverhang 100",
        extra=lambda w, input: f"--sjdbGTFfile {input.annotation} --sjdbOverhang 100",
    log:
        f"logs/genome/star_index.{REF_BUILD}.log",
    cache: True
    wrapper:
        f"{WRAPPER_PATH}/bio/star/index"


rule bowtie2_index:
    input:
        ref=f"resources/{REF_BUILD}/genome.fasta",
    output:
        multiext(
            f"resources/{REF_BUILD}/genome",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        f"logs/genome/{REF_BUILD}/bowtie2_index.{REF_BUILD}.log",
    threads: 8
    cache: True
    wrapper:
        f"{WRAPPER_PATH}/bio/bowtie2/build"


rule hisat2_index_extra:
    input:
        gtf=f"resources/{REF_BUILD}/genome.gtf",
    output:
        multiext(
            f"resources/{REF_BUILD}/genome", ".exons.bed", ".splice_sites.bed"
        ),
    log:
        f"logs/genome/{REF_BUILD}/hisat2_index.features.log",
    conda:
        # relative path to this snakemake file
        workflow.source_path("../envs/hisat2.yaml")
    threads: 8
    # cache: True
    shell:
        r"""
        hisat2_extract_exons.py {input.gtf} > {output[0]}
        hisat2_extract_splice_sites.py {input.gtf} > {output[1]}
        """


rule hisat2_index:
    input:
        fasta=f"resources/{REF_BUILD}/genome.fasta",
        exons=f"resources/{REF_BUILD}/genome.exons.bed",
        splice_sites=f"resources/{REF_BUILD}/genome.splice_sites.bed",
    output:
        multiext(
            f"resources/{REF_BUILD}/genome",
            ".1.ht2",
            ".2.ht2",
            ".3.ht2",
            ".4.ht2",
            ".5.ht2",
            ".6.ht2",
            ".7.ht2",
            ".8.ht2",
        ),
    params:
        prefix=f"resources/{REF_BUILD}/genome",
        extra=lambda w, input: f"--exon {input.exons} --ss {input.splice_sites}",
    log:
        f"logs/genome/{REF_BUILD}/hisat2_index.log",
    threads: 8
    # cache: True
    wrapper:
        f"{WRAPPER_PATH}/bio/hisat2/index"


rule kallisto_index:
    input:
        fasta=f"resources/{REF_BUILD}/genome.cdna.fasta",
    output:
        index=f"resources/{REF_BUILD}/kallisto_genome/transcriptome.idx",
    params:
        extra="",
    log:
        f"logs/genome/{REF_BUILD}/kallisto_index.{REF_BUILD}.log",
    threads: 8
    cache: True
    wrapper:
        f"{WRAPPER_PATH}/bio/kallisto/index"


rule salmon_decoys:
    input:
        genome=f"resources/{REF_BUILD}/genome.fasta",
        transcriptome=f"resources/{REF_BUILD}/genome.cdna.fasta",
    output:
        # Forced, named output files are required in wrapper, cache ignored
        gentrome=f"resources/{REF_BUILD}/genome.transcriptome.fasta",
        decoys=f"resources/{REF_BUILD}/genome.decoys.txt",
    threads: 4
    log:
        f"logs/genome/salmon_decoys.{REF_BUILD}.log",
    cache: False
    wrapper:
        f"{WRAPPER_PATH}/bio/salmon/decoys"


rule salmon_index:
    input:
        sequences=f"resources/{REF_BUILD}/genome.transcriptome.fasta",
        decoys=f"resources/{REF_BUILD}/genome.decoys.txt",
    output:
        directory(f"resources/{REF_BUILD}/salmon_genome"),
    threads: 8
    log:
        f"logs/genome/{REF_BUILD}/salmon_index.{REF_BUILD}.log",
    cache: True
    wrapper:
        f"{WRAPPER_PATH}/bio/salmon/index"
