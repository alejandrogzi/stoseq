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

    TRIMMOMATIC(fastqs, out).set { trimmed_fastqs }

    if (!file("${out}/star/SAindex").exists()) {
        genome_dir = INDEX(fasta, gtf).genome_dir
    } else {
        genome_dir = "${out}/star"
    }

    STAR(fasta, trimmed_fastqs, genome_dir, out).set { star_out }
    // STAR2PASS(fasta, reads_ch, STAR.out.junctions.collect(), genome_dir)
}