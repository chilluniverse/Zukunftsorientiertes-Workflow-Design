#!/usr/bin/env nextflow

//? path definitions
params.genome = "$baseDir/data/references/RNAseq/dmel/*.fasta"    // transcriptome reference files
params.reads = "$baseDir/data/fasta/RNAseq/**/*.fasta"            // fasta raw sequences
params.outdir = "$baseDir/results/RNAseq"                                  // output directory


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!            PROCESSES            !!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//> Build Bowtie2 Index
process BOWTIE2_BUILD {
    tag "$fasta"
    label 'process_high'

    input:
    path(fasta)

    output:
    path('bowtie2')

    script:
    """
    mkdir bowtie2
    bowtie2-build --threads ${params.max_cpus} $fasta bowtie2/index
    """
}

//> Align Sequence
process BOWTIE2_ALIGN {
    maxForks 1
    memory '30 GB'
    errorStrategy 'retry'
    maxErrors 5

    tag "$fasta_file"

    input:
    path fasta_file
    path index_dir

    output:
    path "${fasta_file.simpleName}.sam", emit: sam

    script:
    """
    bowtie2 --very-sensitive-local -N 1 --threads ${params.max_cpus} -f -x ${index_dir}/index -U ${fasta_file} -S "${fasta_file.simpleName}.sam"
    """
}

//> Convert SAM to BAM and Index
process PROCESS_SAM {
    maxForks 1
    memory '30 GB'
    errorStrategy 'retry'
    maxErrors 5

    tag "$sam_file"
    publishDir "$params.outdir/BAM", mode: 'copy', overwrite: 'true', pattern: "{*.tsv,*.stats}"

    input:
    path sam_file

    output:
    path("${sam_file.simpleName}.sorted.bam")
    path("${sam_file.simpleName}.sorted.bam.bai")
    path("${sam_file.simpleName}.sorted.bam.tsv")   ,emit: read_count
    path("${sam_file.simpleName}.bam.stats")

    """
    mkdir samtools

    samtools view       -@ ${params.max_cpus} -Sb ${sam_file.simpleName}.sam -o ${sam_file.simpleName}.bam

    samtools sort       -@ ${params.max_cpus} -m 1300M --output-fmt bam --output-fmt-option nthreads=${params.max_cpus} -o ${sam_file.simpleName}.sorted.bam "${sam_file.simpleName}.bam"

    samtools index      -@ ${params.max_cpus} ${sam_file.simpleName}.sorted.bam -o ${sam_file.simpleName}.sorted.bam.bai

    samtools idxstats   -@ ${params.max_cpus} ${sam_file.simpleName}.sorted.bam &> ${sam_file.simpleName}.sorted.bam.tsv

    samtools stats      -@ ${params.max_cpus} ${sam_file.simpleName}.bam &> ${sam_file.simpleName}.bam.stats
    """
}

//> Build Matrix from read counts
process BUILD_MATRIX {

    publishDir "$params.outdir/DEG", mode: 'copy', overwrite: 'true'

    input:
    val files

    output:
    path "countData.tsv"

    script:
    """
    build_matrix.py --files "$files" --name "countData"
    """
}

//> Generate Metadata file
process GENERATE_METADATA {
    publishDir "$params.outdir/DEG", mode: 'copy', overwrite: 'true'

    input:
    path sam_files

    output:
    path 'metadata.csv'

    script:
    """
    echo "Run,Organism,Time_point" > metadata.csv
    for file in ${sam_files}; do
        filename=\$(basename "\$file" .sorted.bam.tsv)
        organism=\$(echo \$filename | cut -d'_' -f1)
        time=\$(echo \$filename | cut -d'_' -f2)
        sample=\$(echo \$filename | cut -d'_' -f3)

        # Map the time to the correct letter
        case \$time in
            72h)
                letter="A"
                ;;
            96h)
                letter="B"
                ;;
            120h)
                letter="C"
                ;;
            *)
                letter="Unknown"
                ;;
        esac

        echo "\${filename},\${organism},\${letter}_\${time}" >> metadata.csv
    done
    """
}

process DEG {
    publishDir "$params.outdir", mode: 'copy', overwrite: 'true'

    input:
    path countData
    path metaData

    output:
    path ('results')
    path ('results/coseq-clusters_all.csv') ,emit: cluster_all

    script:
    """
    mkdir -p results{plots,upregulated_genes,significantGenes,clustering}
    RNA-seq.r $countData $metaData
    """
}

process EXTRACT_CLUSTER {
    publishDir "$params.outdir", mode: 'copy', overwrite: 'true'

    input:
    path cluster_data

    output:
    path('cluster') ,emit: cluster
    path 'cluster.job'

    script:
    """
    mkdir cluster
    coseq_go.py --clusters $cluster_data
    """
}

//!!!!!!!!!!!!!!!!!!!!!!!!
//!     WORKFLOW
//!!!!!!!!!!!!!!!!!!!!!!!!

workflow {
    //> Alignment / Mapping
    if (params.align) {
        genome = file(params.genome)                    // get transcriptome file from path
        bowtie_index = BOWTIE2_BUILD(genome)            // build bowtie index from transcriptome
        reads = Channel.fromPath(params.reads)          // get reads from path
        BOWTIE2_ALIGN(reads, bowtie_index).sam          // Alignment
                                        | collect
                                        | flatten
                                        | collate( 1 )
                                        | set {bowtie_output}

        read_count = PROCESS_SAM(bowtie_output).read_count 
                                               | collect(flat: false) // Process SAM Files
    } else {
        // if alignment is skipped, get read count data from results-path
        read_count = Channel.fromPath("${params.outdir}/BAM/*.sorted.bam.tsv").collect(flat: false) 
    }

    //> DEG Analysis
    if (params.DEG) {   
        countData = BUILD_MATRIX(read_count)                    // Build Count Matrix for DEG-Analysis
        metaData = GENERATE_METADATA(read_count)                // Build Metadata File for for DEG-Analysis
        cluster_data = DEG(countData, metaData).cluster_all     // Execute DEG Analysis
        EXTRACT_CLUSTER(cluster_data)                           // Write each sub cluster from DEG/coseq-clustering-analsys to unique file and generate .job file for metascape
    }
}