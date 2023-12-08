process QC_RNA {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::muon=0.1.5"
    container "heylf/muon:0.1.5"

    input:
    tuple val(meta), val(samples), val(h5s), val(demuxs)

    output:
    tuple val(meta), path("*.html"), emit: htmls
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: 
    def sample_names = samples.join(",")
    def h5s = h5s.join(",")
    def demux = demuxs.join(",")

    """
    qc_rna.py \\
        --samples $sample_names \\
        --files $h5s

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        qc_rna: \$(qc_rna.py --version | sed 's/qc_rna.py //g')
    END_VERSIONS
    """
}