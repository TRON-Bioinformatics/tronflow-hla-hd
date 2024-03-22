process BAM2FASTQ {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy", pattern: "*.txt"
    module params.bowtie2_module

    conda (params.enable_conda ? "bioconda::samtools=1.18 bioconda::gatk4=4.2.5.0" : null)

    input:
    tuple val(name), val(bam)
    val(contigs)

    output:
    tuple val("${name}"), path("${name}.hla.1.fastq.gz"), path("${name}.hla.2.fastq.gz")

    script:
    """
    # gets reads in the provided regions
    # only non duplicated reads
    samtools view \
        -b \
        -L ${contigs} \
        -F 1024 \
        ${bam} > ${name}.mhc.bam

    # Extract unmap reads
    samtools view -b -f 4 $bam > ${name}.unmap.bam

    #Merge bam files
    samtools merge -o ${name}.merge.bam ${name}.unmap.bam ${name}.mhc.bam

    #Create fastq
    gatk SamToFastq --VALIDATION_STRINGENCY SILENT -I ${name}.merge.bam -F ${name}.hlatmp.1.fastq -F2 ${name}.hlatmp.2.fastq

    #Change fastq ID
    cat ${name}.hlatmp.1.fastq |awk '{if(NR%4 == 1){O=\$0;gsub("/1"," 1",O);print O}else{print \$0}}' | gzip > ${name}.hla.1.fastq.gz
    cat ${name}.hlatmp.2.fastq |awk '{if(NR%4 == 1){O=\$0;gsub("/2"," 2",O);print O}else{print \$0}}' | gzip > ${name}.hla.2.fastq.gz
    """
}
