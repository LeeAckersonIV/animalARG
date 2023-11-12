#!/usr/bin/env nextflow

/* 
 * 'animal-ARG' - A Nextflow pipeline for estimating antimicrobial resistance gene abundance
 */



/*
 * The following pipeline parameters specify the reference genomes
 * and read pairs and can be provided as command line options
 */

params.SRA_IDs = "USER_Defined_List"
params.CARD_db = "$baseDir/data/DB/CARD_NR.99.fasta"
params.GG_db = "$baseDir/data/DB/greenGenes_NR.99.fasta.gz"
params.geneLengthsARG = "$baseDir/data/DB/amr.gene.length.99"
params.geneLengths16S = "$baseDir/data/DB/geneLengths_16srRNA_GreenGenes.99"
params.outdir = "$baseDir/results"
params.QC_files = "$baseDir/QCmetrics"
params.reads = "$baseDir/data/reads/*_{1,2}.fastq.gz"
params.clean_reads = "$baseDir/data/cleanReads/"
params.metadata = "$baseDir/metadata"
params.readDir = "$baseDir/data/reads/"


log.info """\
     A N I M A L  -  A R G     
===============================
SRA IDs      : $params.SRA_IDs
Reads        : $params.reads
ARG Database : $params.CARD_db
16S Database : $params.GG_db
Metadata     : $params.metadata
QC Analysis  : $params.QC_files
Results      : $params.outdir
"""

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

 /* 
 * Import modules 
 */
include { 
  FASTQC;
  fastq_qualityControl;
  FASTQC_on_clean;
  map_ARG; 
  map_16S;
  count_ARG;
  count_16S} from './modules.nf' 

/* 
 * main bioinformatic pipeline 
 */

workflow {

        read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true ) /* this ouputs in tuple format of: [sampleID, [file1, file2]] */
        read_pairs_ch.view()

        // Part 1: Quality Control
        fastqc_ch = FASTQC( read_pairs_ch )
        fastqc_ch.view()
        qc_ch = fastq_qualityControl( read_pairs_ch )
        qc_ch.view()
        fastqc_clean_ch = FASTQC_on_clean( qc_ch )
        fastqc_clean_ch.view()


        // Part 2: Gene Mapping
        bam_ARG_ch = map_ARG( qc_ch ) /* dev-script4.nf */
        bam_ARG_ch.view()
        bam_16S_ch = map_16S( qc_ch ) /* dev-script5.nf */
        bam_16S_ch.view()

        // Part 3: Gene Quantification
        quant_ARG_ch = count_ARG ( bam_ARG_ch )
        quant_ARG_ch.view()
        quant_16S_ch = count_16S ( bam_16S_ch )
        quant_16S_ch.view()

}

