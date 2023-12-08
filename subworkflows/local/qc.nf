//
// Check input samplesheet and get read channels
//

include { QC_RNA } from '../../modules/local/qcrna/main'

workflow QC {
    take:
        h5s

    main:
        ch_versions = Channel.empty()

        QC_RNA( h5s )
        ch_versions = ch_versions.mix(QC_RNA.out.versions)

    emit:
        ch_versions
        htmls = QC_RNA.out.htmls

}