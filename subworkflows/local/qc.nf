//
// Check input samplesheet and get read channels
//

include { QC_RNA } from '../../modules/local/qcrna/main'
include { QC_ATAC } from '../../modules/local/qcatac/main'

workflow QC {
    take:
        h5s

    main:
        ch_versions = Channel.empty()

        QC_RNA( h5s )
        ch_versions = ch_versions.mix(QC_RNA.out.versions)

        QC_ATAC( h5s )
        ch_versions = ch_versions.mix(QC_ATAC.out.versions)

    emit:
        ch_versions
        rna_htmls = QC_RNA.out.htmls
        atac_htmls = QC_ATAC.out.htmls

}
