/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQ_PREPROCESSING                } from '../../../subworkflows/local/fastq_preprocessing/'
include { METAGENOME_ASSEMBLING              } from '../../../subworkflows/local/metagenome_assembling/'
include { POSTPROCESSING_METAGENOME_ASSEMBLY } from '../../../subworkflows/local/postprocessing_metagenome_assembly/'
include { QC_METAGENOME_ASSEMBLY             } from '../../../subworkflows/local/qc_metagenome_assembly/'
include { MULTIQC_SUMMARY                    } from '../../../subworkflows/local/multiqc_summary/'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MG_ASSEMBLER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    
    //
    // PREPROCESSING
    //
    FASTQ_PREPROCESSING(ch_samplesheet)

    // формируем канал нужного формата
    ch_prepared_reads = FASTQ_PREPROCESSING.out.prepared_fastqs.map {meta, illumina_reads, illumina_reads_unpaired ->
        [meta, illumina_reads, illumina_reads_unpaired, [], []]
        }

    //
    // METAGENOME_ASSEMBLING
    //
    METAGENOME_ASSEMBLING(ch_prepared_reads)

    //
    // POSTPROCESSING
    //
    POSTPROCESSING_METAGENOME_ASSEMBLY(METAGENOME_ASSEMBLING.out.contigs)

    //
    // FINAL_QC
    //
    QC_METAGENOME_ASSEMBLY(POSTPROCESSING_METAGENOME_ASSEMBLY.out.prepared_contigs)

    //
    // Collecting QC & versions data
    //
    ch_multiqc_files = ch_multiqc_files.mix(
        FASTQ_PREPROCESSING.out.ch_multiqc_files,
        QC_METAGENOME_ASSEMBLY.out.ch_multiqc_files
    )

    ch_versions = ch_versions.mix(
        FASTQ_PREPROCESSING.out.ch_versions,
        METAGENOME_ASSEMBLING.out.ch_versions,
        POSTPROCESSING_METAGENOME_ASSEMBLY.out.ch_versions,
        QC_METAGENOME_ASSEMBLY.out.ch_versions
    )

    //
    // FINAL_MULTIQC_REPORT
    //
    MULTIQC_SUMMARY(
                    ch_multiqc_files,
                    ch_versions
                   )


    emit:
    
    multiqc_report = MULTIQC_SUMMARY.out.multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
