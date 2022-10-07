#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


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


process HLAHD {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy"

    input:
    tuple val(name), val(fastq1), val(fastq2)

    output:
    tuple val("${name}"), path("*final*")

    script:
    """
    mkdir temp

    # HLA-HD wants its own binaries and bowtie2 in path
    export PATH=${params.hlahd_folder}/bin/:${params.bowtie2_folder}:$PATH

    # HLA HD does not accept gzipped fastq files, unzip them first
    hlahd.sh \
        -m ${params.read_length} \
        -t ${task.cpus} \
        -f ${params.hlahd_folder}/freq_data/ \
        <(zcat ${fastq1}) \
        <(zcat ${fastq2}) \
        ${params.hlahd_folder}/HLA_gene.split.txt \
        ${params.hlahd_folder}/dictionary/ \
        ${name} \
        temp

    # moves the final result to base folder
    mv temp/**/*final* .

    # deletes temp folder
    rm -rf temp
    """
}

workflow {
    HLAHD(input_files)
}
