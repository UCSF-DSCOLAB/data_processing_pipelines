process GTF_GENE_FILTER {
    tag "$fasta"

    input:
    path fasta
    path gtf

    output:
    path "*.gtf"       , emit: genes_gtf

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    filter_gtf_genes.py \\
        --gtf $gtf \\
        --fasta $fasta \\
        --output ${fasta.baseName}_genes.gtf
    """
}