#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BAM2FASTQ } from './modules/00_bam2fastq'
include { HLAHD } from './modules/01_hlahd'


def helpMessage() {
    log.info params.help_message
}

if (params.help) {
    helpMessage()
    exit 0
}

// checks required inputs
if (params.input_fastqs) {
  Channel
    .fromPath(params.input_fastqs)
    .splitCsv(header: ['name', 'fastq1', 'fastq2'], sep: "\t")
    .map{ row-> tuple(row.name, file(row.fastq1), file(row.fastq2)) }
    .set { input_fastqs }
} else if (params.input_bams) {
  Channel
    .fromPath(params.input_bams)
    .splitCsv(header: ['name', 'bam'], sep: "\t")
    .map{ row-> tuple(row.name, file(row.bam)) }
    .set { input_bams }
} else {
  exit 1, "Provide either --input_fastqs or --input_bams"
}


if (params.reference == 'hg38') {
    contigs = "$baseDir/references/contigs_hla_reads_hg38.bed"
}
else if (params.reference == 'hg19') {
    contigs = "$baseDir/references/contigs_hla_reads_hg19.bed"
}
else {
    log.error "--reference only supports hg38 or hg19"
    exit 1
}



workflow {

    if (input_bams) {
        BAM2FASTQ(input_bams, contigs)
        input_fastqs =BAM2FASTQ.out
    }

    HLAHD(input_fastqs)
}
