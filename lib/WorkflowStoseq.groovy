//
// This file holds the WorkflowStoseq class, which is specific to stosseq.nf pipeline
//

class WorkflowStoseq {

    // Citation string for pipeline
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name}, for your analysis please cite:\n\n" +
        "${workflow.manifest.citation}"
    }

    // Generate help string
    public static String help(workflow, params) {
        return "Usage: nextflow run ${workflow.manifest.name} [options]\n\n" +
        "Options:\n" +
        "  --runlist       List of SRA/ENA sample accessions. One per line.\n" +
        "  --genome        Genome to align reads to. Can be specified multiple times.\n" +
        "  --dir           Input directory (if user wants to align their own samples). Default: ${params.dir}\n" +
        "  --out           Output directory. Default: ${params.dir}/results\n" +
        "  --help          Print this message.\n" +
        "  --version       Print version information.\n" +
        "\n" +
        "Examples:\n\n" +
        " If you want to run the pipeline on a list of SRA/ENA accessions, create a file with one accession per line and run:\n\n" +
        "   nextflow run ${workflow.manifest.name} --runlist runlist.txt --genome hg38.fa --gtf hg38.gtf" +
        "\n\n" +
        " If you want to run the pipeline on your own fastq files, place them in a directory and run:\n\n" + 
        "   nextflow run ${workflow.manifest.name} --dir /path/to/fastq --genome hg38.fa --gtf hg38.gtf" +
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
    }
}