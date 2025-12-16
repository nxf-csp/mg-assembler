include { MERGE_LANES } from '../../../modules/local/merge_lanes'
include { FASTQC as FASTQC_RAW } from '../../../modules/nf-core/fastqc'
include { FASTQC as FASTQC_TRIMMED } from '../../../modules/nf-core/fastqc'
include { FASTQC as FASTQC_TRIMMED_UNPAIRED } from '../../../modules/nf-core/fastqc'
include { FASTP } from '../../../modules/local/fastp' 

workflow FASTQ_PREPROCESSING {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_versions
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    // Объединяем данные с разных дорожек
    MERGE_LANES(ch_samplesheet)

    // Первичный QC FASTQ
    FASTQC_RAW(MERGE_LANES.out.fastqs)
    
    // Тримминг, фильтрация, дедупликация
    FASTP(MERGE_LANES.out.fastqs)

    // Вторичный QC
    FASTQC_TRIMMED(FASTP.out.reads)
    FASTQC_TRIMMED_UNPAIRED(FASTP.out.reads_unpaired)

    // Формирование данных о версиях и списка QC файлов
    ch_multiqc_files = ch_multiqc_files.mix(
        FASTQC_RAW.out.zip.collect{it[1]},
        FASTP.out.json.collect{it[1]},
        FASTQC_TRIMMED.out.zip.collect{it[1]},
        FASTQC_TRIMMED_UNPAIRED.out.zip.collect{it[1]}
        )

    ch_versions = ch_versions.mix(
        FASTQC_RAW.out.versions.first(),
        FASTP.out.versions
        )

    emit:
    ch_versions
    ch_multiqc_files
    prepared_fastqs  = FASTP.out.reads.join(FASTP.out.reads_unpaired, by: 0)
}