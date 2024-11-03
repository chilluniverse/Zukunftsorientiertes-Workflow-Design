#!/usr/bin/env nextflow

//? path definitions
params.genome = "$baseDir/data/references/RNAseq/*.fasta"    // transcriptome reference files
params.reads = "$baseDir/data/fasta/RNAseq/**/*.fasta"       // fasta raw sequences
params.outdir = "$baseDir/data/results/RNAseq"               // output directory

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!            PROCESSES            !!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//> Build Bowtie2 Index
//* in preparation for Alignment with Bowtie2
process BOWTIE2_BUILD {
    container 'staphb/bowtie2:2.5.4'

    tag "$fasta"

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
//* using Bowtie2
process BOWTIE2_ALIGN {
    container 'staphb/bowtie2:2.5.4'

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
//* using SAMtools; additionally generate Alignment Statistics from Input BAM
process PROCESS_SAM {
    container 'staphb/samtools:1.20'

    errorStrategy 'retry'
    maxErrors 5

    tag "$sam_file"

    publishDir "$params.outdir/1_alignedReads", mode: 'move', overwrite: 'true', pattern: "{*.stats}"
    publishDir "$params.outdir/1_alignedReads/read_counts", mode: 'copy', overwrite: 'true', pattern: "{*.tsv}"

    input:
    path sam_file

    output:
    path("${sam_file.simpleName}.sorted.bam.tsv")   ,emit: read_count
    path("${sam_file.simpleName}.bam.stats")

    """
    mkdir samtools
    samtools view       -@ ${params.max_cpus} -Sb ${sam_file.simpleName}.sam -o ${sam_file.simpleName}.bam
    samtools sort       -@ ${params.max_cpus} -m 1300M --output-fmt bam --output-fmt-option nthreads=${params.max_cpus} -o ${sam_file.simpleName}.sorted.bam "${sam_file.simpleName}.bam"
    samtools index      -@ ${params.max_cpus} ${sam_file.simpleName}.sorted.bam -o ${sam_file.simpleName}.sorted.bam.bai
    samtools idxstats   -@ ${params.max_cpus} -X ${sam_file.simpleName}.sorted.bam ${sam_file.simpleName}.sorted.bam.bai &> ${sam_file.simpleName}.sorted.bam.tsv
    samtools stats      -@ ${params.max_cpus} ${sam_file.simpleName}.bam &> ${sam_file.simpleName}.bam.stats
    """
}

//> Build Matrix from read counts
//* in preparation for DEG with DESeq2
process BUILD_MATRIX {
    container 'rnaseq:1.0'

    publishDir "$params.outdir/1_Matrix_Metadata_for_DEG", mode: 'copy', overwrite: 'true'

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
//* in preparation for DEG with DESeq2
process GENERATE_METADATA {
    publishDir "$params.outdir/1_Matrix_Metadata_for_DEG", mode: 'copy', overwrite: 'true'

    input:
    path file

    output:
    path 'metadata.csv'

    script:
    """
    echo "Run,Organism,Time_point" > metadata.csv
    for file in ${file}; do
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
//> differential expression gene analysis and clustering
//* with DESeq2 and cosiq
process DEG {
    container 'rnaseq:1.0'

    publishDir "$params.outdir", mode: 'copy', overwrite: 'true'

    input:
    path countData
    path metaData

    output:
    path ("2_DEG")
    path ("2_DEG/clustering/coseq-clusters_all.csv") ,emit: cluster_all

    script:
    """
    mkdir -p 2_DEG/{plots,upregulated_genes,significantGenes,clustering}
    RNA-seq.r $countData $metaData 2_DEG
    """
}

//> Split Genes according to their cluster
process EXTRACT_CLUSTER {
    container 'rnaseq:1.0'

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

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!            WORKFLOW            !!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
workflow {
    //> Alignment / Mapping
    //? execute only if parameter --align is set via CLI on execution
    if (params.align) {
        genome = file(params.genome)                    // get transcriptome file from path
        bowtie_index = BOWTIE2_BUILD(genome)            // build bowtie index from transcriptome
        reads = Channel.fromPath(params.reads)          // get reads from path
        bowtie_output = BOWTIE2_ALIGN(reads, bowtie_index).sam // Alignment
                                        | collect              // collect all Outputs in a list
                                        | flatten              // separate each entry
                                        | collate( 1 )         // publish each entry seperatly to following process

        read_count = PROCESS_SAM(bowtie_output).read_count      // Process SAM Files
                                               | collect(flat: false) // collect all Outputs in a list
    } else {// if alignment is skipped, get read count data from results-path
        read_count = Channel.fromPath("${params.outdir}/0_alignedReads/**/*.sorted.bam.tsv").collect(flat: false) 
    }

    //> DEG Analysis
    //? execute only if parameter --DEG is set via CLI on execution
    if (params.DEG) {   
        countData = BUILD_MATRIX(read_count)                    // Build Count Matrix for DEG-Analysis
        metaData = GENERATE_METADATA(read_count)                // Build Metadata File for for DEG-Analysis
        cluster_data = DEG(countData, metaData).cluster_all     // Execute DEG Analysis
        EXTRACT_CLUSTER(cluster_data)                           // Write each sub cluster from DEG/coseq-clustering-analsys to unique file and generate .job file for metascape
    }
}