//
// This file holds the WorkflowStoseq class, which is specific to stosseq.nf pipeline
//

class WorkflowMain {

    // Citation string for pipeline
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name}, for your analysis please cite:\n\n" +
        "${workflow.manifest.citation}"
    }

    // Just print a dashed line
    public static String dashedLine() {
        return "     ---------------------------------------     "
    }

    // Generate help string
    public static String help(workflow, params) {
        return "Usage: nextflow run ${workflow.manifest.name} [args] [options]\n\n" +
        "Arguments:\n" +
        "  --runlist       List of SRA/ENA sample accessions. One per line.\n" +
        "  --dir           Input directory (if user wants to align their own samples). Default: ${params.dir}\n\n" +
        "Options:\n" +
        "  --genome        Genome to align reads to. Optional if an index already exists in ${params.out}\n" +
        "  --gtf           GTF file to build a genome index with a gene model. Optional if an index already exists in ${params.out}\n" +
        "  --transcriptome Transcriptome to align reads to. Only required if running pipeline with salmon on\n" +
        "  --out           Output directory. Default: ${params.dir}/results\n" +
        "  --help          Print this message.\n" +
        "  --version       Print version information.\n\n" +
        "Flags:\n" +
        "  --fqc           FastQC flag. If set to true, it will produce a quality assesment over trimmed fastqs. Default: false\n" +
        "  --mqc           MultiQC flag. If set to true, it will produce a quality assesment over trimmed fastqs and STAR stats. Default: false\n" +
        "\n" +
        "Examples:\n\n" +
        " If you want to run the basic pipeline on a list of SRA/ENA accessions, create a file with one accession per line and run:\n\n" +
        "   nextflow run ${workflow.manifest.name} --runlist runlist.txt --genome hg38.fa --gtf hg38.gtf" +
        "\n\n" +
        " If you want to run the basic pipeline on your own fastq files, place them in a directory and run:\n\n" + 
        "   nextflow run ${workflow.manifest.name} --dir /path/to/fastq --genome hg38.fa --gtf hg38.gtf" +
        "\n\n" +
        " If you want to run the basic pipeline to quantify transcripts:\n\n" + 
        "   nextflow run ${workflow.manifest.name} --runlist runlist.txt --genome hg38.fa --gtf hg38.gtf --transcriptome hg38_transcriptome.fa --salmon true" +
        "\n\n" +
        "For any bug/question/features please visit the GitHub repository: ${workflow.manifest.homePage}\n\n" +
        "\n\n"
    }

    // Stoseq runner
    public static void run(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(workflow, params)
            System.exit(0)
        }
        

        // Print workflow version and exit on --version
        if (params.version) {
            log.info "${workflow.manifest.name} ${workflow.manifest.version}"
            System.exit(0)
        }

        // Check that either --runlist or --dir is specified
        if (!params.runlist && !params.dir) {
            log.info "\n"
            log.error "Either --runlist or --dir must be specified"
            log.info help(workflow, params)
            System.exit(1)
        }

        def idx = new File("${params.out}/index/SAindex")
        if (!idx.exists() && !params.genome && !params.gtf) {
            log.info "\n"
            log.error "Genome index not found. Please specify --genome and --gtf"
            log.info "\n\n"
            log.info help(workflow, params)
            System.exit(1)
        }

        // if (params.salmon && !params.genome && !params.gtf) {
        //     log.info "\n"
        //     log.error "Genome index not found. Please specify --genome and --gtf"
        //     log.info "\n\n"
        //     log.info help(workflow, params)
        //     System.exit(1)
        // }

        // if (params.salmon && params.genome && params.gtf && !params.transcriptome) {
        //     log.info "\n"
        //     log.error "Transcriptome not found. Please specify --transcriptome"
        //     log.info "\n\n"
        //     log.info help(workflow, params)
        //     System.exit(1)
        // }
    }

    // Stoseq end point
    public static void end(workflow, params, log) {

        if (workflow.success) {
            log.info "\n" 
            log.info dashedLine()
            log.info "\nPipeline finished successfully!" 
            log.info "\n"
            log.info "Completed at: ${workflow.complete}" 
            log.info "Duration    : ${workflow.duration}" 
            log.info "Success     : ${workflow.success}" 
            log.info "workDir     : ${params.out}" 
            log.info "Exit status : ${workflow.exitStatus}"
        } else {
            log.info "\n" 
            log.info dashedLine()
            log.info "\nFailed: ${workflow.errorReport}" 
            log.info "Exit status: ${workflow.exitStatus}".stripIndent()
        }
    }
}
