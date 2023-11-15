#!/usr/env/bin nextflow

process SALMON {

    publishDir "${out}/salmon/${sample}", mode: 'copy', overwrite: 'false'
    cpus 10

    input:
        file transcriptome
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
    -a ${bam} \\
    --threads ${task.cpus} \\
    -numErrorBins '6' --incompatPrior '0.0' \\
    --biasSpeedSamp '5' --fldMax '1000' --fldMean '250' \\
    --fldSD '25' --forgettingFactor '0.65' \\
    --maxReadOcc '100' --numBiasSamples '2000000' \\
    --numAuxModelSamples '5000000' --numPreAuxModelSamples '5000'  \\
    --numGibbsSamples '0'  --numBootstraps '0'  \\
    --thinningFactor '16'  --sigDigits '3' --vbPrior '1e-05' \\
    -o salmon_output

    find ../../ -name ${bam} -type f -delete
    """
}
