process QC_ATAC {
    tag "$meta.id"
    label 'process_medium'

    container "heylf/snapatac2:2.5.1"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "QC_ATAC module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

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
    def demuxs = demuxs.join(",")
    def demuxsArg = demuxs ? "--demux $demuxs" : ""
    def tmpdir = params.tmpdir_atac ? "--tmpdir ${params.tmpdir_atac}" : ""

    """   
    qc_atac.py \\
        --samples $sample_names \\
        --files $h5s \\
        --genome hg38 \\
        $tmpdir \\
        $demuxsArg \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        qc_rna: \$(qc_atac.py --version | sed 's/qc_atac.py //g')
    END_VERSIONS
    """
}