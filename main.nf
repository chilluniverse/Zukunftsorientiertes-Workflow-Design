#!/usr/bin/env nextflow

//? input paths
// fasta raw sequences
reads_basepath = "$baseDir/data/1_fasta/RNAseq"
params.reads = ""
params.dmel_reads = "$reads_basepath/dmel/*.fasta"
params.dmau_reads = "$reads_basepath/dmau/*.fasta"

// transcriptome reference files
references_basepath = "$baseDir/data/0_references/RNAseq" 
params.dmel_genome = "$references_basepath/dmel/*.fasta"
params.dmau_genome = "$references_basepath/dmau/*.fasta"

//? output paths
// output dir
params.outdir = "$baseDir/data/results"

//bowtie2 index
params.dmel_index_prefix = "$baseDir/data/0_index/RNAseq/dmel"
params.dmau_index_prefix = "$baseDir/data/0_index/RNAseq/dmau"

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

process BOWTIE2_ALIGN {
    tag "$fasta_file"

    input:
    path fasta_file
    path(index_dir)

    output:
    path("aligned_${fasta_file.simpleName}.sam")

    script:
    """
    bowtie2 --very-sensitive-local -N 1 --threads ${params.max_cpus} -f -x ${index_dir}/index -U ${fasta_file} -S aligned_${fasta_file.simpleName}.sam
    """
}

workflow {
    if (!params.skipAlignment)
    {
        dmel_genome = file(params.dmel_genome)
        bowtie_index = BOWTIE2_BUILD(dmel_genome)
        bowtie_index.view()
        dmel_reads = Channel.fromPath(params.dmel_reads)
        BOWTIE2_ALIGN(dmel_reads, bowtie_index)
    }
}
