#!/usr/env/bin nextflow

include { QUANT2PASS } from "../modules/quant2pass"
include { SALMON } from "../modules/salmon"


workflow QUANT {
    take:
        fastqs 
        junctions
        genome_dir
        out
        transcriptome

    main:
        QUANT2PASS(fastqs, junctions, genome_dir, out)
        SALMON(transcriptome, QUANT2PASS.out.quant2transcriptome, out)

    emit:
        SALMON.out.quant
}