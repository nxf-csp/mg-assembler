/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    МОДУЛЬ: CONTIG_POSTPROCESSING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Постпроцессинг контигов: приведение имени к новому формату, фильтрация по длине, добавление метаданных
    Входы: FASTA-файл (в т.ч. в .gz архиве)
    Выходы: FASTA-файл с отфильтрованными контигами; FASTA-файл с контигами, не прошедшими фильтрацию
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
process CONTIG_POSTPROCESSING {
    tag "$meta.id"
    label 'process_low'
    stageInMode 'symlink'

    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/biopython:1.84"

    input:
    tuple val(meta), path(contigs)
    val(min_len)
    val(contig_basename)
    val(description)
    val(additional_meta)

    output:
    tuple val(meta), path('*.processed.fa.gz'), emit: processed_contigs
    path  "versions.yml"                                 , emit: versions

    script:
    """
    # Распаковка gzip файлов
    gunzip -f *.fa*.gz 2>/dev/null || true

    # Обработка
    process_contigs.py ${min_len} ${contig_basename} ${description} ${additional_meta}

    # Запаковка
    gzip *.processed.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}