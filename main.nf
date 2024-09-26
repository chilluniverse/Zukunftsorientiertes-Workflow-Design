#!/usr/bin/env nextflow
// ?
//? input paths
// transcriptome reference files
references_basepath = "$baseDir/data/0_references/RNAseq" 
params.dmel_genome = "$references_basepath/dmel/*.fasta"
params.dmau_genome = "$references_basepath/dmau/*.fasta"

// fasta raw sequences
reads_basepath = "$baseDir/data/1_fasta/RNAseq"
params.dmel_reads = "$reads_basepath/dmel/*.fasta"
params.dmau_reads = "$reads_basepath/dmau/*.fasta"
params.reads = "$reads_basepath/**/*.fasta"

//? output paths
// output dir
params.outdir = "$baseDir/data/results_mapping"


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

    samtools sort       -@ ${params.max_cpus} -m 1300M --output-fmt bam --output-fmt-option nthreads=${params.max_cpus} -o ${sam_file.simpleName}.sorted.bam "${sam_file.simpleName}.bam"

    samtools index      -@ ${params.max_cpus} ${sam_file.simpleName}.sorted.bam -o ${sam_file.simpleName}.sorted.bam.bai

    samtools idxstats   -@ ${params.max_cpus} ${sam_file.simpleName}.sorted.bam &> ${sam_file.simpleName}.sorted.bam.tsv

    samtools stats      -@ ${params.max_cpus} ${sam_file.simpleName}.bam &> ${sam_file.simpleName}.bam.stats
    """
}

//> Build Matrix from read counts
process BUILD_MATRIX {

    publishDir "$params.outdir", mode: 'copy', overwrite: 'true'

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
    publishDir "$params.outdir", mode: 'copy', overwrite: 'true'

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
    path ('DEG')

    script:
    """
    mkdir DEG
    RNA-seq.r $countData $metaData
    """
}

//!!!!!!!!!!!!!!!!!!!!!!!!
//!     WORKFLOW
//!!!!!!!!!!!!!!!!!!!!!!!!

workflow {
    if (params.align)
    {
        if (params.dmel)
        {
            genome = params.dmel_genome
            reads = params.dmel_reads
        }
        else if(params.dmau)
        {
            genome = params.dmau_genome
            reads = params.dmau_reads
        }
        else
        {
            genome = params.dmel_genome
            reads = params.reads
        }
        genome = file(genome)
        bowtie_index = BOWTIE2_BUILD(genome)

        reads = Channel.fromPath(reads)
        BOWTIE2_ALIGN(reads, bowtie_index).sam
                                        | collect
                                        | flatten
                                        | collate( 1 )
                                        | set {bowtie_output}

        read_count = PROCESS_SAM(bowtie_output).read_count 
                                               | collect(flat: false)
    } 
    else 
    {
        read_count = Channel.fromPath("${params.outdir}/*.sorted.bam.tsv").collect(flat: false)
    }

    if (params.DEG)
    {   
        countData = BUILD_MATRIX(read_count)

        metaData = GENERATE_METADATA(read_count)

        DEG(countData, metaData)
    }
    
}
