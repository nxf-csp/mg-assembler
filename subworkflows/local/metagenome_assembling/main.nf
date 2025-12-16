include { SPADES } from '../../../modules/local/spades'

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

    ch_versions.mix(SPADES.out.versions)

    emit:
    contigs  = SPADES.out.contigs
    versions = ch_versions

}