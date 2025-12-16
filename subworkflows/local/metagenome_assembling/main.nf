include { SPADES } from '../../../modules/local/spades'
include { CONTIG_POSTPROCESSING } from '../../../modules/local/contig_postprocessing'

workflow METAGENOME_ASSEMBLING {
    take:
    ch_reads // канал вида [мета, illumina_reads, illumina_reads_unpaired, pacbio_reads, nanopore_reads]

    main:
    // Дефолтные значения
    yml          = []
    hmm          = []
    isolate      = false
    metagenome   = true
    single_cell  = false
    rna          = false
    rna_viral    = false
    plasmid      = false
    meta_plasmid = false
    meta_viral   = false
    bio          = false
    corona       = false
    sewage       = false

    ch_versions = channel.empty()

    SPADES(
           ch_reads,
           yml,
           hmm,
           isolate,
           metagenome,
           single_cell,
           rna,
           rna_viral,
           plasmid,
           meta_plasmid,
           meta_viral,
           bio,
           corona,
           sewage
          )

    CONTIG_POSTPROCESSING(
                          SPADES.out.contigs,
                          params.ctg_min_len,
                          params.ctg_basename,
                          params.ctg_description,
                          params.ctg_additional_meta,
                         )

    ch_versions.mix(
        SPADES.out.versions,
        CONTIG_POSTPROCESSING.out.versions
    )

    emit:
    prepared_contigs = CONTIG_POSTPROCESSING.out.processed_contigs
    versions         = ch_versions

}