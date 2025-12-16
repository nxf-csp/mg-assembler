include { CONTIG_POSTPROCESSING } from '../../../modules/local/contig_postprocessing'

workflow POSTPROCESSING_METAGENOME_ASSEMBLY {
    take:
    ch_contigs // канал вида [мета, fastas]

    main:
    ch_versions = channel.empty()
    
    CONTIG_POSTPROCESSING(
                          ch_contigs,
                          params.ctg_min_len,
                          params.ctg_basename,
                          params.ctg_description,
                          params.ctg_additional_meta,
                         )

    ch_versions.mix(CONTIG_POSTPROCESSING.out.versions)

    emit:
    prepared_contigs = CONTIG_POSTPROCESSING.out.processed_contigs
    versions         = ch_versions
}