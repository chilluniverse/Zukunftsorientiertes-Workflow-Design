#!/usr/bin/env nextflow

//? path definitions
params.download_file = "$baseDir/data/fasta/SRA_Download.csv"
params.data_path = "$baseDir/data"

params.rna_ref = "https://ftp.flybase.org/genomes/dmel/dmel_r6.54_FB2023_05/fasta/dmel-all-transcript-r6.54.fasta.gz"
params.atac_ref = "https://ftp.flybase.org/genomes/dmel/dmel_r6.54_FB2023_05/fasta/dmel-all-chromosome-r6.54.fasta.gz"
params.atac_gtf = "http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.54_FB2023_05/gtf/dmel-all-r6.54.gtf.gz"
params.chip = ["http://furlonglab.embl.de/labData/publications/2012/Junion-et-al-2012_Cell/signal_files/Pnr_4-6h_allIPvsallM_log2ratio_5probes-smoothed_Junion2012.bigwig", "http://furlonglab.embl.de/labData/publications/2012/Junion-et-al-2012_Cell/signal_files/Pnr_6-8h_allIPvsallM_log2ratio_5probes-smoothed_Junion2012.bigwig"]
params.docker = "$baseDir/docker/*/"

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!            PROCESSES            !!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//> Download Read from SRA
//* and save with simple file name
process READS {
    conda 'bioconda::sra-tools=3.1.1 conda-forge::ossuuid=1.6.2'

    tag "$rename"

    publishDir "$params.data_path/fasta/$path", mode: 'move', overwrite: 'true'

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

//> Download ChIP-Chip data
//* and convert to BED file format
process CHIP_CHIP {
    conda 'bioconda::ucsc-bigwigtobedgraph'

    publishDir "$params.data_path/fasta/ChIP-chip", mode: 'move', overwrite: 'true'

    input:
    val(link)

    output:
    path('*.bed')

    shell:
    """
    wget $link
    link=$link
    filename=\${link##*/}
    bigWigToBedGraph \$filename \${filename%.bigwig}.bed
    """
}

//> Build Docker Container
//* with dependencies that are not available as containers 
process BUILD_CONTAINER {
    input:
    path folder

    output:
    stdout

    script:
    """
    name=$folder
    docker buildx build -t \${name,,}:1.0 $folder
    """   
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!            WORKFLOW            !!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

workflow {

    CHIP_CHIP(Channel.fromList(params.chip))

    REFERENCE_RNA(params.rna_ref)

    REFERENCE_ATAC(params.atac_ref, params.atac_gtf)

    download_file = Channel
                        .fromPath(params.download_file)
                        .splitCsv(header: true)

    download_file.map{ row -> tuple(row.run, row.rename, row.path) } | READS

    docker = Channel.fromPath(params.docker, type: 'dir')

    params.noDocker = false
    if ( !params.noDocker ) {
        BUILD_CONTAINER(docker) | view
    }
}