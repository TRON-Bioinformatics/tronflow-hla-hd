process BAM2FASTQ {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy", pattern: "*.txt"
    module params.bowtie2_module

    conda (params.enable_conda ? "bioconda::samtools=1.18" : null)

    input:
    tuple val(name), val(bam)
    val(contigs)

    output:
    tuple val("${name}"), path("${name}.hla.1.fastq.gz"), path("${name}.hla.2.fastq.gz")

    script:
    """
    samtools index $bam

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
    rm -f ${name}.unmap.bam
    rm -f ${name}.mhc.bam

    # sort/collate BAM by read name and convert to FASTQ
     samtools collate -u -O ${name}.merge.bam | samtools fastq \
        -1 ${name}.hlatmp.1.fastq -2 ${name}.hlatmp.2.fastq -s /dev/null -0 /dev/null

    #Change fastq ID
    awk '{if(NR%4 == 1){O=\$0;gsub("/1"," 1",O);print O}else{print \$0}}' ${name}.hlatmp.1.fastq | gzip > ${name}.hla.1.fastq.gz
    awk '{if(NR%4 == 1){O=\$0;gsub("/2"," 2",O);print O}else{print \$0}}' ${name}.hlatmp.2.fastq | gzip > ${name}.hla.2.fastq.gz

    rm -f ${name}.hlatmp.1.fastq
    rm -f ${name}.hlatmp.2.fastq
    """
}
