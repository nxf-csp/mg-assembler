include { MULTIQC                } from '../../../modules/nf-core/multiqc/main'
include { paramsSummaryMultiqc   } from '../../nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../../nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../local/utils_nfcore_mg-assembler_pipeline'
include { paramsSummaryMap       } from 'plugin/nf-schema'

workflow MULTIQC_SUMMARY {
    take:
    ch_multiqc_files
    ch_versions

    main:

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'mg-assembler_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
    ch_multiqc_files.collect(),
    ch_multiqc_config.toList(),
    ch_multiqc_custom_config.toList(),
    ch_multiqc_logo.toList(),
    [],
    []
)

    emit:
    
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
}