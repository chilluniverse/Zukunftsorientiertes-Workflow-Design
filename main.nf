#!/usr/bin/env nextflow
// ?
//? input paths
// transcriptome reference files
references_basepath = "$baseDir/data/0_references/RNAseq" 
params.dmel_genome = "$references_basepath/dmel/*.fasta"
params.dmau_genome = "$references_basepath/dmau/*.fasta"

// fasta raw sequences
reads_basepath = "$baseDir/data/1_fasta/RNAseq"
params.reads = ""
params.dmel_reads = "$reads_basepath/dmel/*.fasta"
params.dmau_reads = "$reads_basepath/dmau/*.fasta"

//? output paths
// output dir
params.outdir = "$baseDir/data/results_mapping"

//bowtie2 index
params.dmel_index_prefix = "$baseDir/data/0_index/RNAseq/dmel"
params.dmau_index_prefix = "$baseDir/data/0_index/RNAseq/dmau"

//!!!!!!!!!!!!!!!!!!!!!!!!
//!     PROCESSES
//!!!!!!!!!!!!!!!!!!!!!!!!

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
    memory '20 GB'

    tag "$fasta_file"

    input:
    path fasta_file
    path index_dir

    output:
    path "${fasta_file.simpleName}.sam"

    script:
    """
    bowtie2 --very-sensitive-local -N 1 --threads ${params.max_cpus} -f -x ${index_dir}/index -U ${fasta_file} -S "${fasta_file.simpleName}.sam"
    """
}

//> Convert SAM to BAM and Index
process PROCESS_SAM {

    tag "$sam_file"
    publishDir "$params.outdir", mode: 'copy', overwrite: 'true'

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

    samtools sort       -@ ${params.max_cpus} -m 1G --output-fmt bam --output-fmt-option nthreads=${params.max_cpus} -o ${sam_file.simpleName}.sorted.bam "${sam_file.simpleName}.bam"

    samtools index      -@ ${params.max_cpus} ${sam_file.simpleName}.sorted.bam -o ${sam_file.simpleName}.sorted.bam.bai

    samtools idxstats   -@ ${params.max_cpus} ${sam_file.simpleName}.sorted.bam &> ${sam_file.simpleName}.sorted.bam.tsv

    samtools stats      -@ ${params.max_cpus} ${sam_file.simpleName}.bam &> ${sam_file.simpleName}.bam.stats
    """
}

//!!!!!!!!!!!!!!!!!!!!!!!!
//!     WORKFLOW
//!!!!!!!!!!!!!!!!!!!!!!!!

workflow {
    if (!params.skipAlign)
    {
        dmel_genome = file(params.dmel_genome)
        bowtie_index = BOWTIE2_BUILD(dmel_genome)

        dmel_reads = Channel.fromPath(params.dmel_reads)
        bowtie_output = BOWTIE2_ALIGN(dmel_reads, bowtie_index)

        read_count = PROCESS_SAM(bowtie_output).read_count 
                                               | collect(flat: false)
                                               | view
    }


}
