process BAM2FASTQ {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}/${name}", mode: "copy", pattern: "*.txt"
    module params.bowtie2_module

    conda (params.enable_conda ? "bioconda::samtools=1.18 bioconda::bedtools=2.31.1" : null)

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

    # sort BAM by read name
    samtools sort -n ${name}.merge.bam -o ${name}.merge.sorted.bam
    rm -f ${name}.merge.bam

    #Create fastq
    bedtools bamtofastq -i ${name}.merge.sorted.bam -fq ${name}.hlatmp.1.fastq -fq2 ${name}.hlatmp.2.fastq

    #Change fastq ID
    cat ${name}.hlatmp.1.fastq |awk '{if(NR%4 == 1){O=\$0;gsub("/1"," 1",O);print O}else{print \$0}}' | gzip > ${name}.hla.1.fastq.gz
    cat ${name}.hlatmp.2.fastq |awk '{if(NR%4 == 1){O=\$0;gsub("/2"," 2",O);print O}else{print \$0}}' | gzip > ${name}.hla.2.fastq.gz

    rm -f ${name}.hlatmp.1.fastq
    rm -f ${name}.hlatmp.2.fastq
    """
}
