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
//* nextflow run main.nf --runlist --genome --gtf --out [to download runs and index]
//*
//* nextflow run main.nf --dir --genome --gtf --out [to use local fastqs and index]


// modules
include { TRIMMOMATIC } from './modules/trim'
include { INDEX } from './modules/index'
include { STAR } from './modules/star'
include { STAR2PASS } from './modules/star2pass'

// subworkflows
include { GET_FASTQ } from './subworkflows/get'


WorkflowStoseq.run(workflow, params, log)


workflow {

    if (params.dir == '.') {
        fastqs = GET_FASTQ(params.runlist).fastqs
    } else {
        fastqs = Channel.fromFilePairs("${params.dir}/*_{1,2}.fastq*.gz")
    }

    out = file(params.out ?: "${params.dir}/results")
    out.mkdir()

    genome_dir = "${out}/index"

    TRIMMOMATIC(fastqs, out, params.dir)

    if (!file("${out}/index/SAindex").exists()) {
        genome = file(params.genome)
        gtf = file(params.gtf)
        INDEX(genome, gtf, out)
    } 
    
    STAR(TRIMMOMATIC.out.trim_paired, genome_dir, out).set { star_out }
    STAR2PASS(TRIMMOMATIC.out.trim_paired, STAR.out.junctions.collect(), genome_dir, out)
}