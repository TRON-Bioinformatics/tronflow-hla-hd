#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HLAHD } from './modules/01_hlahd'


def helpMessage() {
    log.info params.help_message
}

if (params.help) {
    helpMessage()
    exit 0
}

// checks required inputs
if (params.input_files) {
  Channel
    .fromPath(params.input_files)
    .splitCsv(header: ['name', 'fastq1', 'fastq2'], sep: "\t")
    .map{ row-> tuple(row.name, file(row.fastq1), file(row.fastq2)) }
    .set { input_files }
} else {
  exit 1, "Input file not specified!"
}


workflow {
    HLAHD(input_files)
}
