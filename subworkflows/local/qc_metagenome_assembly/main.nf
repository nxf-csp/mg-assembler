include { SEQKIT_STATS } from '../../../modules/nf-core/seqkit/stats'

workflow QC_METAGENOME_ASSEMBLY {
    take:
    ch_contigs // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    SEQKIT_STATS(ch_contigs)

    // Формирование данных о версиях и списка QC файлов
    ch_multiqc_files = ch_multiqc_files.mix(SEQKIT_STATS.out.stats)

    ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions)

    emit:
    versions  = ch_versions
    mqc_files = ch_multiqc_files

}