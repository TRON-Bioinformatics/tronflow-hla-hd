process HLAHD {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy", pattern: "*.txt"
    module params.bowtie2_module

    input:
    tuple val(name), val(fastq1), val(fastq2)

    output:
    tuple val("${name}"), path("*final*")

    script:
    """
    mkdir temp

    # HLA-HD wants its own binaries and bowtie2 in path
    export PATH=${params.hlahd_folder}/bin/:${params.bowtie2_folder}:\$PATH
    export LD_LIBRARY_PATH=${params.ld_library_path}

    zcat ${fastq1} > input_fastq1.fastq
    zcat ${fastq2} > input_fastq2.fastq

    # HLA HD does not accept gzipped fastq files, unzip them first
    hlahd.sh \
        -m ${params.read_length} \
        -t ${task.cpus} \
        -f ${params.hlahd_folder}/freq_data/ \
        input_fastq1.fastq \
        input_fastq2.fastq \
        ${params.hlahd_folder}/HLA_gene.split.txt \
        ${params.hlahd_folder}/dictionary/ \
        ${name} \
        temp

    # moves the final result to base folder
    mv temp/${name}/result/* .

    # deletes temp folder
    rm -rf temp
    rm -f input_fastq1.fastq
    rm -f input_fastq2.fastq
    """
}
