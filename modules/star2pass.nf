#!/usr/env/bin nextflow


process STAR2PASS {
    
    publishDir "${out}/star-2pass/${sample}", mode: 'copy', overwrite: 'false'
    cpus 15

    input:
        tuple val(sample), file(fastqs)
        path junctions
        path genome_dir
        val out

    output:
        tuple val(sample), file("*.final.out"), emit: bam2stats
        tuple val(sample), file("*ReadsPerGene*.tab"), emit: counts
        // tuple val(sample), file("*sortedByCoord*.bam"), emit: bam2pass
        // tuple val(sample), file("*.bai"), emit: bam2index
        // tuple val(sample), file("*toTranscriptome*.bam"), emit: bam2transcriptome

    script:
    """
    mkdir -p ${out}/star-2pass/${sample}

    STAR --genomeDir ${genome_dir} \\
    --readFilesCommand zcat \\
    --readFilesIn ${fastqs} \\
    --runThreadN ${task.cpus} \\
    --outFileNamePrefix "${sample}.2pass." \\
    --sjdbFileChrStartEnd ${junctions} \\
    --quantMode GeneCounts TranscriptomeSAM \\
    --outSAMtype BAM SortedByCoordinate \\
    --outSAMattributes NH HI AS nM MD ch \\
    --outSAMattrRGline ID:${sample} LB:library PL:illumina PU:machine SM:${sample} \\
    --outSAMunmapped None

    # samtools index "${sample}.2pass.Aligned.sortedByCoord.out.bam"

    rm *.bam
    """
}