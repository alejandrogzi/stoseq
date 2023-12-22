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
include { MAKE_TRANSCRIPTOME } from './modules/to-trans/to-trans'
include { TRIMMOMATIC } from './modules/trim/trimmomatic'
include { TRIM_GALORE } from './modules/trim/trim_galore'
include { RRNA } from './modules/sortmerna/rrna'
include { FASTQC } from './modules/fastqc/fastqc'
include { INDEX } from './modules/star/index'
include { STAR } from './modules/star/star'
include { STAR2PASS } from './modules/star/star2pass'

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
    gentrome_dir = "${out}/salmon"

    if (params.trim) {
      trim_fq = TRIM_GALORE(fastqs, params.out, params.dir)
    } else {
      trim_fq = TRIMMOMATIC(fastqs, params.out, params.dir)
    }

    rrna_db = Channel.from("${params.rrna}/*.fasta").map{fa -> file(fa)}.collect()
    rrna_fq = RRNA(trim_fq.trim_paired, rrna_db)

    if (params.fqc) {
        FASTQC(TRIMMOMATIC.out.trim_paired, params.out)
    }

    if (!file("${params.out}/salmon/quant.sf").exists()) {
        genome = file(params.genome)
        gtf = file(params.gtf)
        transcriptome = MAKE_TRANSCRIPTOME(genome, gtf)
        salmon_idx = SALMON_INDEX(genome, transcriptome, params.out)
    }

    if (!file("${params.out}/index/SAindex").exists()) {
        idx = INDEX(genome, gtf, params.out)
    } 

    bam = STAR(rrna_fq.reads, genome_dir, out)
    bam2pass = STAR2PASS(rrna_fq.reads, bam.junctions.collect(), genome_dir, out)
    salmon = QUANT(bam2pass.bam, gentrome_dir, out)

}

workflow.onComplete{ 
 WorkflowMain.end(workflow, params, log) 
}
