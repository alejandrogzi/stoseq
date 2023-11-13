#!/usr/env/bin nextflow


include { TRIMMOMATIC } from './modules/trim'
include { INDEX } from './modules/index'
include { STAR } from './modules/star'
include { STAR2PASS } from './modules/star2pass'


workflow {
    params.dir = 'NONE'
    params.out = null

    if (params.dir == 'NONE') {
        println "Please provide a directory with the --dir flag"
        System.exit(1)
    }

    out = file(params.out ?: "${params.dir}/results")
    out.mkdir()

    fastqs = Channel.fromFilePairs("${params.dir}/*_{1,2}.fastq*.gz")
    fasta = file("${params.dir}/*.fa")
    gtf = file("${params.dir}/*.gtf")
    genome_dir = "${out}/index"

    TRIMMOMATIC(fastqs, out).set { trimmed_fastqs }

    if (!file("${out}/index/SAindex").exists()) {
        INDEX(fasta, gtf, out)
    } 
    
    STAR(fasta, trimmed_fastqs, genome_dir, out).set { star_out }
    STAR2PASS(fasta, trimmed_fastqs, STAR.out.junctions.collect(), genome_dir, out)
}