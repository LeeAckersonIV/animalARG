/*
 * Process 1A: Run FastQC
 */
process FASTQC {
    tag "${sample_id}" 
    publishDir params.QC_files

    input:
    tuple val(sample_id), file(reads)

    output:
    path "${sample_id}_fastqc.log"

    script:
    """
    fastqc --outdir $params.QC_files ${reads[0]} ${reads[1]} > ${sample_id}_fastqc.log
    """
 }

/*
 * Process 1B: Quality Control ~ Adapter Trimming, Quality Filtering, Read Trimming
 */
 process fastq_qualityControl {
    tag "QC on $sample_id"
    publishDir params.clean_reads

        input:
    tuple val(sample_id), file(reads)


        output:
    tuple val("$sample_id"), path("${sample_id}_*.clean_*.fastq.gz")

    script:
    """
    bbmap/bbduk.sh in1=${reads[0]} in2=${reads[1]} out1=${sample_id}_1.clean_1.fastq.gz out2=${sample_id}_2.clean_2.fastq.gz outm=${sample_id}_1.outm_1.fastq.gz outm2=${sample_id}_2.outm_2.fastq.gz outs=${sample_id}.outs.fastq.gz ref=bbmap/resources/adapters.fa qtrim=r trimq=20 minlen=50
    """
}

/*
 * Process 1C: Run FastQC on Clean Reads
 */
process FASTQC_on_clean {
    tag "${sample_id}"
    publishDir params.QC_files

    input:
    tuple val(sample_id), file(reads)

    output:
    path "${sample_id}_clean_fastqc.log"

    script:
    """
    fastqc --outdir $params.QC_files ${reads[0]} ${reads[1]} > ${sample_id}_clean_fastqc.log
    """
 }

/*
 * Process 2A: Align to ARG Reference: CARD
 */
process map_ARG {
    tag "Align ${sample_id} -> ARG Database"
    publishDir params.outdir
    input:
    tuple val(sample_id), path(clean_reads)

        output:
    tuple val("$sample_id"), path("${sample_id}_ARG.bam")


    script:
    """
    bwa mem -t 8 ${params.CARD_db} ${clean_reads[0]} ${clean_reads[1]} 2> ${sample_id}_ARG_bwaMem.err | samtools view -bS - > ${sample_id}_ARG.bam 2> ${sample_id}_ARG_samtoolsView.err
    """
}

/*
 * Process 2B: Align to 16S rRNA Reference: GreenGenes 
 */
process map_16S {
    tag "Align ${sample_id} -> 16S Database"
    publishDir params.outdir
    input:
    tuple val(sample_id), path(clean_reads)

    output:
    tuple val("$sample_id"), path("${sample_id}_16S.bam")

    script:
    """
    bwa mem -t 8 ${params.GG_db} ${clean_reads[0]} ${clean_reads[1]} 2> ${sample_id}_16S_bwaMem.err | samtools view -bS - > ${sample_id}_16S.bam 2> ${sample_id}_16S_samtoolsView.err
    """
}

 /*
 * Process 3A: Gene Quantification ~ ARGs
 */
process count_ARG {
	tag "Quantify ARG Abundance: ${sample_id}"
	publishDir params.outdir
        
    input: 
    tuple val(sample_id), path(bam)
        
    output: 
    tuple val("$sample_id"), path("${sample_id}_absAbun_ARG")

    shell:
    '''
    join -t $'\t' -a 1 -a 2 -e - -o '0,1.2,2.2,1.3,2.3' <(samtools view -q 1 -f 0x40 !{bam} | cut -f 1,3,5 | sort -k 1,1) <(samtools view -q 1 -F 0x40 !{bam} | cut -f 1,3,5 | sort -k 1,1) | perl -wne 'chomp $_; @line = split /\t/, $_; if ($line[1] eq $line[2]) { print $line[0], "\t", $line[1], "\n"; } elsif ( $line[1] eq "-" ) { print $line[0], "\t", $line[2], "\n"; } elsif ( $line[2] eq "-" ) { print $line[0], "\t", $line[1], "\n"; } elsif ($line[3] >= $line[4]) { print $line[0], "\t", $line[1], "\n"; } else { print $line[0], "\t", $line[2], "\n"; }' | awk '{print $1"::::"$2}' | sort -k1,1 | bedtools groupby -g 1 -c 1 -o collapse,count | sed 's/::::/\t/' | sort -k1,1 | bedtools groupby -g 1 -c 2,4,4 -o collapse,collapse,sum | perl -wne 'chomp $_; @line = split /\t/, $_; @genes = split /,/, $line[1]; @counts = split /,/, $line[2]; for (my $i = 0; $i <= $#genes; $i++) { print $genes[$i], "\t", $counts[$i]/$line[3], "\n"; }' | sort -k1,1 | bedtools groupby -g 1 -c 2 -o sum > !{sample_id}_absAbun_ARG
    '''
}

/*
 * Process 3B: Gene Quantification ~ 16S rRNA
 */
process count_16S {
	tag "Quantify 16S rRNA Abundance: ${sample_id}"
	publishDir params.outdir
        
    input: 
    tuple val(sample_id), path(bam)
        
    output: 
    tuple val("$sample_id"), path("${sample_id}_absAbun_16S")

    shell:
    '''
    join -t $'\t' -a 1 -a 2 -e - -o '0,1.2,2.2,1.3,2.3' <(samtools view -q 1 -f 0x40 !{bam} | cut -f 1,3,5 | sort -k 1,1) <(samtools view -q 1 -F 0x40 !{bam} | cut -f 1,3,5 | sort -k 1,1) | perl -wne 'chomp $_; @line = split /\t/, $_; if ($line[1] eq $line[2]) { print $line[0], "\t", $line[1], "\n"; } elsif ( $line[1] eq "-" ) { print $line[0], "\t", $line[2], "\n"; } elsif ( $line[2] eq "-" ) { print $line[0], "\t", $line[1], "\n"; } elsif ($line[3] >= $line[4]) { print $line[0], "\t", $line[1], "\n"; } else { print $line[0], "\t", $line[2], "\n"; }' | awk '{print $1"::::"$2}' | sort -k1,1 | bedtools groupby -g 1 -c 1 -o collapse,count | sed 's/::::/\t/' | sort -k1,1 | bedtools groupby -g 1 -c 2,4,4 -o collapse,collapse,sum | perl -wne 'chomp $_; @line = split /\t/, $_; @genes = split /,/, $line[1]; @counts = split /,/, $line[2]; for (my $i = 0; $i <= $#genes; $i++) { print $genes[$i], "\t", $counts[$i]/$line[3], "\n"; }' | sort -k1,1 | bedtools groupby -g 1 -c 2 -o sum > !{sample_id}_absAbun_16S
    '''
}
