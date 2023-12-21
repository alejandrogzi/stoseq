#!/usr/env/bin nextflow

//* stoseq: A portable storage-optimized RNA-seq processing pipeline
//* 
//* author: Alejandro Gonzales-Irribarren
//* email: jose.gonzalesdezavala1@unmsm.edu.pe
//* github: alejandrogzi
//* version: 1.0.0
//*
//*
//* usage:
//*
//* nextflow run main.nf --runlist --genome --gtf --out [to download runs]
//*
//* nextflow run main.nf --dir --genome --gtf --out [to use local fastqs]


// modules
include { TRIMMOMATIC } from './modules/trim'
include { RRNA } from './modules/rrna'
include { FASTQC } from './modules/fastqc'
include { INDEX } from './modules/index'
include { STAR } from './modules/star'
include { STAR2PASS } from './modules/star2pass'

// subworkflows
include { GET_FASTQ } from './subworkflows/get'
include { QUANT } from './subworkflows/salmon-quant'


WorkflowMain.run(workflow, params, log)


workflow {

    if (params.dir == '.') {
        fastqs = GET_FASTQ(params.runlist).fastqs
    } else {
        fastqs = Channel.fromFilePairs("${params.dir}/*_{1,2}.fastq.gz")
    }

    out = file(params.out)
    out.mkdir()
    genome_dir = "${out}/index"

    trim_fq = TRIMMOMATIC(fastqs, params.out, params.dir)

    rrna_db = Channel.from("${params.rrna}/*.fasta").map{fa -> file(fa)}.collect()
    rrna_fq = RRNA(trim_fq.trim_paired, rrna_db)

    if (params.fqc) {
        FASTQC(TRIMMOMATIC.out.trim_paired, params.out)
    }

    if (!file("${params.out}/index/SAindex").exists()) {
        genome = file(params.genome)
        gtf = file(params.gtf)
        idx = INDEX(genome, gtf, params.out)
    } 
    
    bam = STAR(rrna_fq.reads, genome_dir, out)

    if (!params.salmon) {
        bam2pass = STAR2PASS(trim_fq.trim_paired, bam.junctions.collect(), genome_dir, out)
    } else {
        QUANT(TRIMMOMATIC.out.trim_paired, STAR.out.junctions.collect(), genome_dir, params.out, params.transcriptome)
    }
}

workflow.onComplete{ 
 WorkflowMain.end(workflow, params, log) 
}
