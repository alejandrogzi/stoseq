#!/usr/env/bin nextflow

process QUANT {

    publishDir "${out}/salmon", mode: 'copy', overwrite: 'false'
    cpus 10

    input:
        path gtf
        path transcriptome
        tuple val(sample), file(bam)
        val out

    output:
        tuple val(sample), path("${sample}/*genes.sf"), emit: gene
        tuple val(sample), path("${sample}/quant.sf"), emit: transcript

    script:
    """
    mkdir -p ${out}/salmon/${sample}

    salmon quant \\
    -t ${transcriptome} \\
    -l A \\
    --geneMap ${gtf} \\
    -a ${bam} \\
    --threads ${task.cpus} \\
    -o ${sample}

    # find ../../ -name ${bam} -type f -delete
    """
}
