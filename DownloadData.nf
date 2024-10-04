#!/usr/bin/env nextflow

//? path definitions
params.download_file = "$baseDir/data/fasta/SRA_Download.csv"
params.data_path = "$baseDir/data"

params.rna_ref = "https://ftp.flybase.org/genomes/dmel/dmel_r6.54_FB2023_05/fasta/dmel-all-transcript-r6.54.fasta.gz"
params.atac_ref = "https://ftp.flybase.org/genomes/dmel/dmel_r6.54_FB2023_05/fasta/dmel-all-chromosome-r6.54.fasta.gz"
params.atac_gtf = "http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.54_FB2023_05/gtf/dmel-all-r6.54.gtf.gz"

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!            PROCESSES            !!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//> Download Read from SRA
//* and save with simple file name
process READS {
    conda 'bioconda::sra-tools'

    tag "$rename"

    publishDir "$params.data_path/$path", mode: 'move', overwrite: 'true'

    input:
    tuple val(run), val(rename), val(path)

    output:
    path "${rename}.fasta"

    shell:
    """
    fasterq-dump $run --outfile ${rename}.fasta --threads ${params.max_cpus} --progress --mem 2048MB --fasta
    """
}

//> Download RNA reference
//* from provided link
process REFERENCE_RNA {
    publishDir "$params.data_path/references/RNAseq", mode: 'move', overwrite: 'true'

    input:
    val(link)

    output:
    path('*.fasta')

    shell:
    """
    wget $link
    link=$link
    gzip -d \${link##*/}
    """
}

//> Download ATAC reference
//* from provided link and discard mitochondrial genome
process REFERENCE_ATAC {
    conda 'bioconda::seqkit'

    publishDir "$params.data_path/references/ATACseq", mode: 'move', overwrite: 'true'

    input:
    val(link_fasta)
    val(link_gtf)

    output:
    path('*_wo_mito.fasta')
    path('*.gtf')

    shell:
    """
    wget $link_fasta
    link=$link_fasta
    archive=\${link##*/}
    gzip -d \$archive


    seqkit grep -v -r -p 'mitochondrion_genome' \${archive%.gz} -o \${archive%.fasta.gz}_wo_mito.fasta

    wget $link_gtf
    link=$link_gtf
    gzip -d \${link##*/}
    """
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!            WORKFLOW            !!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

workflow {
    REFERENCE_RNA(params.rna_ref)

    REFERENCE_ATAC(params.atac_ref, params.atac_gtf)

    download_file = Channel
                        .fromPath(params.download_file)
                        .splitCsv(header: true)

    download_file.map{ row -> tuple(row.run, row.rename, row.path) } | READS
}