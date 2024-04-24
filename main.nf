#!/usr/bin/env nextflow

//Enable DSL2 syntax
nextflow.enable.dsl = 2

// Define params and set default paths
params.reads = "$baseDir/raw_data/*_{1,2}.fastq.gz"
params.trimmed_dir = "$baseDir/trimmed_output"
params.asm_dir = "$baseDir/asm_output"
params.qa_dir = "$baseDir/qa_output"
params.gntyp_dir = "$baseDir/gntyp_output"

// Check if trimmed output directory exists, if not create it
if (!file(params.trimmed_dir).exists()) {
    file(params.trimmed_dir).mkdirs()
}

// Define input Channel for raw reads
reads_ch = Channel.fromPath( params.reads )

// Process to run BBDUK
process bbduk {
    // Set a tag
    tag "Trim ${sample_id}"

    // Publish the output to trimmed dir
    publishDir params.trimmed_dir, mode:'copy'

    // Input Reads
    input:
    tuple val(sample_id), path(reads) 

    // Output files
    output:
    tuple val(sample_id), path("${sample_id}_trimmed_{1,2}.fastq.gz")

    // Command to run BBDuk
    script:
    """
    bbduk.sh \\
        in=${reads[0]} \\
        in2=${reads[1]} \\
        out=${sample_id}_trimmed_1.fastq.gz \\
        out2=${sample_id}_trimmed_2.fastq.gz \\
        ref=$baseDir/phix_data/sequence.fasta \\
        k=31 \\
        hdist=1 \\
        t=4 \\
        qtrim=r \\
        trimq=20 \\
        minlen=36
    """
}

// Check if assembly output directory exists, if not create it
if (!file(params.asm_dir).exists()) {
    file(params.asm_dir).mkdirs()
}

// Process to run SKESA for assembly
process skesa {
    // Set a tag for the process
    tag "Assemble $sample_id"

    // Publish the output to assembly directory
    publishDir params.asm_dir, mode:'copy'

    // Input trimmed reads
    input:
    tuple val(sample_id), path(trimmed_reads) 
    
    // Output assembly files
    output:
    tuple val(sample_id), path("${sample_id}.fna"), emit: assembled_reads

    // Command to run SKESA
    script:
    """
    skesa --memory 5 --reads ${trimmed_reads[0]} ${trimmed_reads[1]} --contigs_out ${sample_id}.fna
    """
}

// Check if qa output directory exists, if not create it
if (!file(params.qa_dir).exists()) {
    file(params.qa_dir).mkdirs()
}

// Quality Assessment using QUAST
process quast {
    // Set a tag for the process
    tag "QA $sample_id"

    // Publish the output to qa directory
    publishDir params.qa_dir, mode:'copy'

    // Input Assembled files
    input:
    tuple val(sample_id), path(assembled_reads) 

    // Output QA files
    output:
    path "${sample_id}"

    // Command to run QUAST
    script:
    """
    mkdir ${sample_id}
    quast ${assembled_reads} -o ${sample_id}
    """
}

// Check if genotyping output directory exists, if not create it
if (!file(params.gntyp_dir).exists()) {
    file(params.gntyp_dir).mkdirs()
}

// Genotyping with MLST
process mlst {
    // Set a tag for the process
    tag "Genotyping $sample_id"

    // Publish the output to genotyping directory
    publishDir params.gntyp_dir, mode:'copy'

    // Input Assembled files
    input:
    tuple val(sample_id), path(assembled_reads) 

    // Output MLST files
    output:
    path("${sample_id}_MLST_Summary.tsv")

    // Command to run QUAST
    script:
    """
    mlst ${assembled_reads} > ${sample_id}_MLST_Summary.tsv
    """
}


// Define the workflow
workflow {
    // Create a Channel from input file pairs
    reads = Channel.fromFilePairs(params.reads, checkIfExists: true)

    // Trim reads using BBDUK
    trimmed_ch = bbduk( reads )

    // Assemble trimmed reads using SKESA
    assembly_ch = skesa( trimmed_ch )

    // QA using QUAST
    quast( assembly_ch )

    // Genotyping using MLST
    mlst ( assembly_ch )
}