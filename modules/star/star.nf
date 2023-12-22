#!/usr/env/bin nextflow


process STAR {

    publishDir "${out}/star/${sample}", mode: 'copy', overwrite: 'false'
    cpus 15

    input:
        tuple val(sample), file(fastqs)
        path genome_dir
        val out

    output:
        // tuple val(sample), file("*.bam"), emit: bam
        tuple val(sample), file("*.final.out"), emit: bam_stats
        path("*.tab"), emit: junctions
        

    script:
    """
    mkdir -p ${out}/star/${sample}

    STAR --genomeDir ${genome_dir} \\
    --outFileNamePrefix "${sample}." \\
    --readFilesIn ${fastqs} \\
    --runThreadN ${task.cpus} \\
    --readFilesCommand zcat --outFilterType BySJout \\
    --outSAMtype BAM Unsorted --alignSJoverhangMin 5 \\
    --alignSJDBoverhangMin 3 --outFilterMismatchNmax 10 \\
    --outSAMunmapped None

    rm ${sample}.Aligned.out.bam
    """ 
}