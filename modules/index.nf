#!/usr/env/bin nextflow


process INDEX {
    cpus 12

    input:
        file fasta
        file gtf
        val out
    output:
        path(star), emit: genome_dir
    script:
    """
    mkdir star

    STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir $out/star \\
    --genomeFastaFiles $fasta \\
    --genomeSAindexNbases 14 \\
    --sjdbGTFfile $gtf \\
    --sjdbOverhang 100
    """
}