#!/usr/bin/env nextflow

//? path definitions
params.genome = "$baseDir/data/references/ATACseq/*.fasta"    // transcriptome reference files
params.gtf = "$baseDir/data/references/ATACseq/*.gtf"
// params.reads = "$baseDir/data/fasta/ATACseq/**/*.fasta"            // fasta raw sequences
params.reads = "$baseDir/data/fasta/tmp/*.fasta"            // fasta raw sequences
params.chip = "$baseDir/data/fasta/ChIP-chip/*.bed"
params.pnr = "$baseDir/data/pnr-motif/*.motif"

params.outdir = "$baseDir/results/ATACseq"                         // output directory

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!            PROCESSES            !!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//> Build Bowtie2 Index
process BOWTIE2_BUILD {
    container 'staphb/bowtie2:2.5.4'
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
    container 'staphb/bowtie2:2.5.4'
    maxForks 1
    memory '30 GB'
    errorStrategy 'retry'
    maxErrors 5

    tag "$fasta_file"

    input:
    path fasta_file
    path index_dir

    output:
    tuple val("${fasta_file.simpleName}"), path('mapping/') , emit: mapping

    script:
    """
    mkdir mapping
    bowtie2 --no-unal -X 2000 --threads ${params.max_cpus} -f -x ${index_dir}/index -U ${fasta_file} -S "mapping/${fasta_file.simpleName}.sam" --met-file mapping/${fasta_file.simpleName}.metrics
    """
}

//> Convert SAM to BAM and Index
process PROCESS_SAM {
    maxForks 1
    memory '30 GB'
    errorStrategy 'retry'
    maxErrors 1

    tag "$sample_name"
    publishDir "$params.outdir/BAM", mode: 'copy', overwrite: 'true', pattern: "{*.bed,*.stats}"

    input:
    tuple val(sample_name), path(mapping)

    output:
    path("${sample_name}.sorted.bam")
    path("${sample_name}.sorted.bam.bai")
    path("${sample_name}.bam.stats")
    path("${sample_name}_rmD.bam")
    path("${sample_name}_rmD.bed") ,emit: bed_file

    """
    samtools view       -@ ${params.max_cpus} -h -S -b -o ${sample_name}.bam "$mapping/${sample_name}.sam"

    samtools sort       -@ ${params.max_cpus} -m 1300M --output-fmt bam -o ${sample_name}.sorted.bam ${sample_name}.bam

    samtools index      -@ ${params.max_cpus} ${sample_name}.sorted.bam -o ${sample_name}.sorted.bam.bai

    samtools stats      -@ ${params.max_cpus} ${sample_name}.bam &> ${sample_name}.bam.stats

    picard MarkDuplicates I=${sample_name}.sorted.bam O=${sample_name}_rmD.bam M="$mapping/${sample_name}.metrics" REMOVE_DUPLICATES=TRUE

    bedtools bamtobed -i "${sample_name}_rmD.bam" > "${sample_name}_rmD.bed"

    """
}

//> PEAK CALLING
process PEAK_CALLING {
    container 'community.wave.seqera.io/library/macs2:2.2.9.1--37875900d3e7f01c'

    publishDir "$params.outdir/Peaks/", mode: 'copy', overwrite: 'true'

    tag "${bed_file.simpleName}"

    input:
    path bed_file

    output:
    tuple val("${bed_file.simpleName}"), path("${bed_file.simpleName}/NA_summits.bed"), emit: summit
    path "${bed_file.simpleName}", emit: macs2

    script:
    """
    mkdir ${bed_file.simpleName}
	macs2 callpeak -t $bed_file -g dm --nomodel --shift -100 --extsize 200 -q 0.01 --bdg --outdir ${bed_file.simpleName}
    """
}

//> Annotation
process ANNOTATE {
    tag "$sample"

    publishDir "$params.outdir/Annotation", mode: 'copy', overwrite: 'true'
    
    input:
    tuple val(sample), path(summits_bed)

    output:
    path "${sample}_annotation.csv"

    script:
    """
    annotatePeaks.pl $summits_bed dm6 -gtf $params.gtf > ${sample}_annotation.csv
    """
}

//> Filtering NarrowPeaks
process FILTER_NARROWPEAKS {
    maxForks 2
    memory '15 GB'

    tag "$macs2"

    input:
    path macs2
    path annotation

    output:
    path "${annotation.simpleName}_filtered.csv"    
    path "${annotation.simpleName}.narrowPeak"
    path "${annotation.simpleName}.bed"         ,emit: filtered_peaks

    // cut -f10 ${annotation} | grep -v 'NA'  > "${annotation.simpleName}_filtered.csv"
    script:
    """
    awk -F'\\t' '\$10 != "NA"' ${annotation} > "${annotation.simpleName}_filtered.csv"
    filter_peaks.py --annotation ${annotation.simpleName}_filtered.csv --input "$macs2/NA_peaks.narrowPeak" --output "${annotation.simpleName}.narrowPeak"
    cut -f 1-6 ${annotation.simpleName}.narrowPeak > ${annotation.simpleName}.bed
    """
}

//> de novo Pnr motif search
process GENERATE_MOTIF {
    maxForks 2
    memory '15 GB'

    tag "$bed_file"

    publishDir "$params.outdir/Pnr-Motif/${bed_file}", mode: 'copy', overwrite: 'true'

    input:
    path bed_file

    output:
    path("motif"), emit: motif

    script:
    """
    mkdir motif
    findMotifsGenome.pl ${bed_file} dm3 motif -size 200 -mask -p ${params.max_cpus}
    """
}

process MERGE_PEAKS {
    container 'pegi3s/bedtools:2.31.0'
    maxForks 1
    memory '30 GB'

    tag "$filtered_peaks_bed"

    input:
    path filtered_peaks_bed

    output:
    path("merged_peaks.bed")

    script:
    """
    cat ${filtered_peaks_bed} > "merged_peaks_tmp.bed"
    sort -k 1,1 -k2,2n "merged_peaks_tmp.bed" > "merged_peaks_sorted.bed"
    bedtools merge -i "merged_peaks_sorted.bed"  > "merged_peaks_woID.bed"
    awk -F'\\t' 'BEGIN {OFS=FS} {if (\$4 == "") \$4 = "peak_" counter++; print}' "merged_peaks_woID.bed" > "merged_peaks.bed"
    """
}

process FIND_MOTIF {
    maxForks 2
    memory '15 GB'

    tag "$motif"

    publishDir "$params.outdir/peaksWithPnr", mode: 'copy', overwrite: 'true'

    input:
    path merged_peaks
    path reference
    path motif

    output:
    path ("peakAnalysis.tsv")
    path ("peaksWithPnr.bed"), emit: bed_file

    script:
    """
    mkdir peakAnalysis
    mkdir preparsed
    findMotifsGenome.pl $merged_peaks \
    $reference \
    "peakAnalysis/" \
    -find $motif \
    -p ${params.max_cpus} \
    -cache 2000 \
    -preparse \
    > "peakAnalysis.tsv" 2> "analysis_log.txt"
    cut -f 1 "peakAnalysis.tsv" | sort | uniq -c | wc -l

    filter_peaks.py --annotation peakAnalysis.tsv --input $merged_peaks --output "peaksWithPnr.bed"
"""
}

process ANNOTATE_PNR_PEAKS {
    maxForks 1
    memory '30 GB'

    publishDir "$params.outdir/pnrPeaksAnnotation", mode: 'copy', overwrite: 'true'
    
    input:
    path(summits_bed)

    output:
    path "pnrPeaks_annotation.csv"          , emit: annotation
    path "pnrPeaks_annotation_onlyGenes.txt", emit: gene_list

    script:
    """
    annotatePeaks.pl $summits_bed dm6 -gtf $params.gtf > pnrPeaks_annotation_tmp.csv
    sort -t\$'\\t' -k12,12 -u pnrPeaks_annotation_tmp.csv > pnrPeaks_annotation.csv;cut -f 12 pnrPeaks_annotation.csv > pnrPeaks_annotation_onlyGenes.txt
    """
}

process GET_SIGNIFICANT_GENES {
    maxForks 1
    memory '30 GB'

    input:
    path gene_list
    path annotation

    output:
    path("significant_genes.txt")

    script:
    """
    significant_genes.py --genes $gene_list --annotation $annotation --outname significant_genes.txt
    """
}

//!!!!!!!!!!!!!!!!!!!!!!!!
//!     WORKFLOW
//!!!!!!!!!!!!!!!!!!!!!!!!

workflow {
    //> Alignment / Mapping
    if (params.align) {
        genome = file(params.genome)                    // get genome file from path
        bowtie_index = BOWTIE2_BUILD(genome)            // build bowtie index from transcriptome
        reads = Channel.fromPath(params.reads)          // get reads from path
        bowtie_output = BOWTIE2_ALIGN(reads, bowtie_index).mapping         // Alignment
                                                            | collect
                                                            | flatten
                                                            | collate( 2 )

        bed_file = PROCESS_SAM(bowtie_output).bed_file 
                                               | collect
                                               | flatten
    } else {
        // if alignment is skipped, get read count data from results-path
        bed_file = Channel.fromPath("${params.outdir}/BAM/*.bed").collect(flat: false) 
    }

    if (params.peakcalling) {
        summit_bed = PEAK_CALLING(bed_file).summit
        macs2 = PEAK_CALLING.out.macs2
        annotation = ANNOTATE(summit_bed)
        filtered_peaks = FILTER_NARROWPEAKS(macs2, annotation).filtered_peaks
                                                                | collect
    }

    if (params.pnrmotif) {
        bed_files = Channel.fromPath(params.chip)
        GENERATE_MOTIF(bed_files)
    }

    if (params.motif) {
        motif = Channel.fromPath(params.pnr)
        merged_peaks = MERGE_PEAKS(filtered_peaks)
        bed_file = FIND_MOTIF(merged_peaks, genome, motif).bed_file
        pnr_annotation = ANNOTATE_PNR_PEAKS(bed_file).annotation
        gene_list = ANNOTATE_PNR_PEAKS.out.gene_list
        GET_SIGNIFICANT_GENES(gene_list, pnr_annotation)
    }

}
