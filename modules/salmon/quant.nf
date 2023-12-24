#!/usr/env/bin nextflow

process QUANT {

    publishDir "${out}/salmon/${sample}", mode: 'copy', overwrite: 'false'
    cpus 10

    input:
        path gtf
        path transcriptome
        tuple val(sample), file(bam)
        val out

    output:
        tuple val(sample), path("*.sf"), emit: quant

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
