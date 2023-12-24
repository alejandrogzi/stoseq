#!/usr/bin/env nextflow

process SALMON_INDEX {

    publishDir "${out}/salmon_index", mode: 'copy'
    cpus 8

    input:
    path genome_fasta
    path transcript_fasta
    val out

    output:
    path "salmon"      , emit: index

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def get_decoy_ids = "grep '^>' $genome_fasta | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 > decoys.txt"
    def gentrome      = "gentrome.fa"
    if (genome_fasta.endsWith('.gz')) {
        get_decoy_ids = "grep '^>' <(gunzip -c $genome_fasta) | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 > decoys.txt"
        gentrome      = "gentrome.fa.gz"
    }
    """
    $get_decoy_ids
    sed -i.bak -e 's/>//g' decoys.txt
    cat $transcript_fasta $genome_fasta > $gentrome

    salmon \\
        index \\
        --threads $task.cpus \\
        -t $gentrome \\
        -d decoys.txt \\
        -i salmon

    """
}
