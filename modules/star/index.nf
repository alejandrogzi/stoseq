#!/usr/env/bin nextflow


process INDEX {
    cpus 12

    input:
        file fasta
        file gtf
        val out

    output:
        path(index), emit: genome_dir
        
    script:
    """
    mkdir index

    STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir $out/index \\
    --genomeFastaFiles $fasta \\
    --genomeSAindexNbases 14 \\
    --sjdbGTFfile $gtf \\
    --sjdbOverhang 100
    """
}