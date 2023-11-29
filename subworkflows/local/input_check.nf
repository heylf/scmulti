//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_h5_channel(it) }
        .set { h5s }

    emit:
    h5s                                       // channel: [ val(meta), [ h5s ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_h5_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample

    // add path(s) of the fastq file(s) to the meta map
    def h5_meta = []
    if (!file(row.h5).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> h5 file does not exist!\n${row.h5}"
    }
    h5_meta = [ meta, file(row.h5) ]
    return h5_meta
}
