#!/usr/env/bin nextflow

// stoseq: A portable storage-optimized RNA-seq processing pipeline
//
// author: Alejandro Gonzales-Irribarren
// email: jose.gonzalesdezavala1@unmsm.edu.pe
// github: alejandrogzi
// version: 1.0.0
//
// Usage:
//
// nextflow run main.nf --runlist --genome --gtf --out [to download runs]
//
// nextflow run main.nf --dir --genome --gtf --out [to use local fastqs]


// modules
include { MAKE_TRANSCRIPTOME } from './modules/to-trans/to-trans'
include { TRIMMOMATIC } from './modules/trim/trimmomatic'
include { TRIM_GALORE } from './modules/trim/trim_galore'
include { RRNA } from './modules/sortmerna/rrna'
include { RIBODETECTOR } from './modules/ribodetector/ribodetector'
include { FASTQC } from './modules/fastqc/fastqc'
include { INDEX } from './modules/star/index'
include { STAR } from './modules/star/star'
include { STAR2PASS } from './modules/star/star2pass'
include { SALMON_INDEX } from './modules/salmon/index'
include { QUANT } from './modules/salmon/quant'

// subworkflows
include { GET_FASTQ } from './subworkflows/get'

WorkflowMain.run(workflow, params, log)


workflow {
    
     dir = params.dir
     if (dir == '.') {
         fastqs = GET_FASTQ(params.runlist).fastqs
     } else {
         fastqs = Channel.fromFilePairs("${dir}/*_{1,2}.f*q.gz")
     }

     out = file(params.out)
     out.mkdir()

     genome = file(params.genome)
     genome_dir = "${out}/index"

     if (params.trim) {
       trim_fq = TRIM_GALORE(fastqs, out)
     } else {
       trim_fq = TRIMMOMATIC(fastqs, out)
     }

     if (params.ribodetector) {
         rfastq = RIBODETECTOR(trim_fq, out)
     } else {
       rrna_db = Channel.from("./assets/rrna_dbs/*.fasta").map{fa -> file(fa)}.collect()
       rfastq = RRNA(trim_fq, rrna_db, out)
     }
 
     if (params.fqc) {
         FASTQC(trim_fq, out)
     }
 
      if (!file("${out}/index/SAindex").exists()) {
         idx = INDEX(genome, gtf, out)
     }

     gtf = file(params.gtf)
     transcriptome = MAKE_TRANSCRIPTOME(genome, gtf)

     bam = STAR(rfastq.reads, genome_dir, out)
     bam2pass = STAR2PASS(rfastq.reads, bam.junctions.collect(), genome_dir, out)
     salmon = QUANT(gtf, transcriptome, bam2pass.bam, out)
}

workflow.onComplete{ 
 WorkflowMain.end(workflow, params, log) 
}
