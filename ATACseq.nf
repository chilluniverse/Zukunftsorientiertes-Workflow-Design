#!/usr/bin/env nextflow

//? path definitions
//@param $baseDir - Path where this Pipeline is located at
params.genome = "$baseDir/data/references/ATACseq/*.fasta"        // genome reference files
params.gtf = "$baseDir/data/references/ATACseq/*.gtf"             // gts reference file
params.reads = "$baseDir/data/fasta/ATACseq/**/*.fasta"           // fasta raw sequences
params.chip = "$baseDir/data/fasta/ChIP-chip/*.bed"               // ChIP-chip bed-files
params.pnr_motif = "$baseDir/data/pnr-motif/*.motif"              // Pnr-motif
params.pnr_motif_opt = "$baseDir/data/pnr-motif/optional/*.motif" // Optional Pnr-motif(s)
params.outdir = "$baseDir/data/results/ATACseq"                   // output directory

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

    label 'half'

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
//* using SAMtools; additionally generate Alignment Statistics from Input BAM
process PROCESS_SAM {
    container 'staphb/samtools:1.20'

    label 'half'

    errorStrategy 'retry'
    maxErrors 1

    tag "$mapping"

    publishDir "$params.outdir/0_alignedReads", mode: 'move', overwrite: 'true', pattern: "{*.stats}"

    input:
    tuple val(sample_name), path(mapping)

    output:
    path("${sample_name}.bam.stats")
    path("${sample_name}.sorted.bam") ,emit: bam_file

    """
    samtools view       -@ ${params.max_cpus} -h -S -b -o ${sample_name}.bam "$mapping/${sample_name}.sam"

    samtools sort       -@ ${params.max_cpus} -m 1300M --output-fmt bam -o ${sample_name}.sorted.bam ${sample_name}.bam

    samtools index      -@ ${params.max_cpus} ${sample_name}.sorted.bam -o ${sample_name}.sorted.bam.bai

    samtools stats      -@ ${params.max_cpus} ${sample_name}.bam &> ${sample_name}.bam.stats
    """
}

//> Remove Duplicate Reads
//* with picard MarkDuplicates, included in gatk
process MARK_DUPLICATES {
    container 'broadinstitute/gatk:4.6.0.0'

    label 'third'

    tag "$in_file"

    publishDir "$params.outdir/0_alignedReads", mode: 'move', overwrite: 'true', pattern: "{*.metrics}"

    input:
    path in_file

    output:
    path("${in_file.simpleName}_rmD.bam"), emit: bam_file
    path("${in_file.simpleName}.metrics")

    script:
    """
    gatk MarkDuplicates -I ${in_file} -O ${in_file.simpleName}_rmD.bam -M "${in_file.simpleName}.metrics" --REMOVE_DUPLICATES TRUE
    """
}

//> Convert BAM in BED File Format
//* using bedtools
process BAM_TO_BED {
    container 'pegi3s/bedtools:2.31.0'

    label 'third'

    tag "$bam_file"

    publishDir "$params.outdir/0_alignedReads/bed", mode: 'copy', overwrite: 'true', pattern: "{*.bed}"

    input:
    path bam_file

    output:
    path ("${bam_file.simpleName}.bed"), emit: bed_file

    script:
    """
    bedtools bamtobed -i "${bam_file}" > "${bam_file.simpleName}.bed"
    """
}

//> PEAK CALLING
//* using macs2
process PEAK_CALLING {
    container 'community.wave.seqera.io/library/macs2:2.2.9.1--37875900d3e7f01c'

    label 'third'

    tag "${bed_file.simpleName}"

    publishDir "$params.outdir/1_Peaks/", mode: 'copy', overwrite: 'true'

    input:
    path bed_file

    output:
    tuple val("${bed_file.simpleName}"), path("${bed_file.simpleName}/NA_summits.bed"), emit: summit
    path "${bed_file.simpleName}/NA_peaks.narrowPeak", emit: narrowPeaks

    script:
    """
    mkdir ${bed_file.simpleName}
	macs2 callpeak -t $bed_file -g dm --nomodel --shift -100 --extsize 200 -q 0.01 --bdg --outdir ${bed_file.simpleName}
    """
}

//> Annotation
//* using annotatePeaks.pl from HOMER
process ANNOTATE {
    container 'atacseq:1.0'

    label 'third'

    tag "$sample"

    publishDir "$params.outdir/2_PeakAnnotation", mode: 'copy', overwrite: 'true'

    input:
    tuple val(sample), path(summits_bed)
    path gtf

    output:
    path "${sample}_annotation.csv"

    script:
    """
    annotatePeaks.pl $summits_bed dm6 -gtf $gtf > ${sample}_annotation.csv
    """
}

//> Filter NarrowPeaks
//* filter out peaks without any valid annotation
process FILTER_NARROWPEAKS {
    container 'atacseq:1.0'

    label 'half'

    tag "$narrowPeaks"

    input:
    path narrowPeaks
    path annotation

    output:
    path "${annotation.simpleName}.bed"         ,emit: filtered_peaks

    script:
    """
    awk -F'\\t' '\$10 != "NA"' ${annotation} > "${annotation.simpleName}_filtered.csv"
    filter_peaks.py --annotation ${annotation.simpleName}_filtered.csv --input "$narrowPeaks" --output "${annotation.simpleName}.narrowPeak"
    cut -f 1-6 ${annotation.simpleName}.narrowPeak > ${annotation.simpleName}.bed
    """
}

//> de novo Pnr motif search
//* Generate de novo Pnr-motif using ChIP-chip data
process GENERATE_MOTIF {
    container 'atacseq:1.0'

    label 'half'

    tag "$bed_file"

    publishDir "$params.outdir/de-novo_Pnr-motif", mode: 'move', overwrite: 'true'

    input:
    path bed_file

    output:
    path("${bed_file.simpleName}"), emit: motif

    script:
    """
    mkdir ${bed_file.simpleName}
    findMotifsGenome.pl ${bed_file} dm3 ${bed_file.simpleName} -size 200 -mask -p ${params.max_cpus}
    """
}

//> Merge all Peak Files
//* into one .bed file
process MERGE_PEAKS {
    container 'pegi3s/bedtools:2.31.0'

    label 'half'

    tag "$filtered_peaks_bed"

    input:
    path filtered_peaks_bed

    output:
    path("merged_peaks.bed")

    script:
    """
    cat ${filtered_peaks_bed} > "merged_peaks_tmp.bed"
    sort -k 1,1 -k2,2n "merged_peaks_tmp.bed" > "merged_peaks_sorted.bed"
    awk -F'\\t' 'BEGIN {OFS=FS} {if (\$4 == "") \$4 = "peak_" counter++; print}' "merged_peaks_sorted.bed" > "merged_peaks.bed"
    """
}

//> Find Pnr-Motif in Peaks
//* using findMotifsGenome.pl from HOMER
process FIND_MOTIF {
    container 'atacseq:1.0'

    label 'half'

    tag "$motif"

    publishDir "$params.outdir/3_Peaks-Pnr", mode: 'copy', overwrite: 'true'

    input:
    path merged_peaks
    path reference
    file motif
    path opt_motif

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
    -opt $opt_motif \
    -p ${params.max_cpus} \
    -cache 2000 \
    -preparse \
    > "peakAnalysis.tsv" 2> "analysis_log.txt"
    cut -f 1 "peakAnalysis.tsv" | sort | uniq -c | wc -l

    filter_peaks.py --annotation peakAnalysis.tsv --input $merged_peaks --output "peaksWithPnr.bed"
    """
}

//> Annotate Peaks containig Pnr-Motif
//* Annotate Peaks to dm6 ref genome
process ANNOTATE_PNR_PEAKS {
    container 'atacseq:1.0'

    label 'half'

    publishDir "$params.outdir/4_Peaks-Pnr-Annotation", mode: 'copy', overwrite: 'true'

    input:
    path(summits_bed)
    path gtf

    output:
    path "pnrPeaks_annotation.csv"          , emit: annotation

    script:
    """
    cat ${summits_bed} > tmp.bed
    annotatePeaks.pl tmp.bed dm6 -gtf $gtf > pnrPeaks_annotation_tmp.csv
    { head -n 1 pnrPeaks_annotation_tmp.csv && tail -n +2 pnrPeaks_annotation_tmp.csv | sort -t\$'\\t' -k12,12; } > pnrPeaks_annotation.csv
    """
}

//> Get significant Genes
//* significant = genes appear certain amount of time (here: 10)
process GET_SIGNIFICANT_GENES {
    container 'atacseq:1.0'

    label 'half'

    publishDir "$params.outdir/4_Peaks-Pnr-Annotation", mode: 'copy', overwrite: 'true', pattern: "{all_genes.txt}"
    publishDir "$params.outdir/5_SignificantGenes", mode: 'move', overwrite: 'true'

    input:
    path annotation
    val count

    output:
    path("all_genes.txt")
    path("significant_genes.txt")

    script:
    """
    cut -f 12 $annotation > all_genes.txt
    sed -i '/^\$/d' all_genes.txt
    significant_genes.py --genes all_genes.txt --annotation $annotation --outname significant_genes.txt --count $count
    """
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!            WORKFLOW            !!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
workflow {
    //> Alignment / Mapping
    genome = file(params.genome)                        // get genome file from path
    //? execute only if parameter --align is set via CLI on execution
    if ( params.align ) {
        bowtie_index = BOWTIE2_BUILD(genome)            // build bowtie index from genome

        reads = Channel.fromPath(params.reads)          // get reads from path
        bowtie_output = BOWTIE2_ALIGN(reads, bowtie_index).mapping          // Alignment
                                                            | collect       // collect all Outputs in a list
                                                            | flatten       // separate each entry
                                                            | collate( 2 )  // publish two consecutive entry at a time to following process

        bam_file = PROCESS_SAM(bowtie_output).bam_file      // Convert SAM to BAM and process
        rmD_bam_file = MARK_DUPLICATES(bam_file).bam_file   // Remove Duplicate Reads
        bed_file = BAM_TO_BED(rmD_bam_file).bed_file        // Convert BAM Files to BED
                                               | collect    // collect all Outputs in a list
                                               | flatten    // publish each item separatly to next process

    } else {
        // if alignment is skipped, get read count data from results-path
        bed_file = Channel.fromPath("${params.outdir}/0_alignedReads/**/*.bed")
    }

    //> Peak Calling
    gtf = file(params.gtf)                              // get gtf reference from path
    //? execute only if parameter --peakcalling is set via CLI on execution
    if ( params.peakcalling ||  params.motif ) {
        summit_bed = PEAK_CALLING(bed_file).summit      // Call Peaks; (1) return bed file
        narrowPeaks = PEAK_CALLING.out.narrowPeaks      // (2) return narrowPeaks file

        annotation = ANNOTATE(summit_bed, gtf)          // Annotate bed files

        filtered_peaks = FILTER_NARROWPEAKS(narrowPeaks, annotation).filtered_peaks // filter out peaks without valid annotation
                                                                | collect   // collect all Outputs in a list
    }

    //> de novo Pnr-Motif
    //? execute only if parameter --pnrmotif is set via CLI on execution
    if ( params.pnrmotif ) {
        bed_files = Channel.fromPath(params.chip)       // get ChIP-Chip files from path
        GENERATE_MOTIF(bed_files)                       // generate de novo Pnr-Motif
    }

    //> Find Pnr-Motif in Peaks
    //? execute only if parameter --motif is set via CLI on execution
    if ( params.motif ) {
        merged_peaks = MERGE_PEAKS(filtered_peaks)                   // merge all Peak Files into one

        pnr_motif = file(params.pnr_motif)                           // get selected Pnr-Motif
        opt_motif = Channel.fromPath(params.pnr_motif_opt).collect() // get more selected Pnr-Motif
        bed_file = FIND_MOTIF(merged_peaks, genome, pnr_motif, opt_motif).bed_file.collect() // find Pnr-Motif in those Peaks

        pnr_annotation = ANNOTATE_PNR_PEAKS(bed_file, gtf).annotation // Annotate Peaks with Pnr-Motif

        num_of_appearances = 10                                   // Define how many reads are mandatory to call a gen significant
        GET_SIGNIFICANT_GENES(pnr_annotation, num_of_appearances) // Get significant genes
    }
}
