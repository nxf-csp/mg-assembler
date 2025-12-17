process MERGE_LANES {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(r1_files), path(r2_files)
    
    output:
    tuple val(meta), path("${meta.id}_*.fastq.gz"), emit: fastqs
    
    script:
    """
    cat ${r1_files.join(' ')} > ${meta.id}_1.fastq.gz
    cat ${r2_files.join(' ')} > ${meta.id}_2.fastq.gz
    """
}