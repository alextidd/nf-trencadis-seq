
// split subsetted BAMs by cell, count coverage
process get_coverage {
  tag "${meta.id}_${meta.gene}"
  label 'week10gb'
  publishDir "${params.out_dir}/runs/${meta.id}/${meta.gene}/",
    mode: "copy",
    pattern: "*_coverage.tsv"

  input:
  tuple val(meta), path(bam), path(bai),
        path(celltypes),
        path(mutations),
        path(cell_barcodes),
        path(gene_bed), path(gene_features), path(exonic_bed)

  output:
  tuple val(meta), path("${meta.id}_${meta.gene}_coverage.tsv"),
        path(celltypes),
        path(cell_barcodes),
        path(gene_bed), path(gene_features), path(exonic_bed),
        emit: cov

  script:
  """
  module load samtools-1.19/python-3.12.0

  # get gene info
  chr=\$(cut -f1 ${gene_bed})
  start=\$(cut -f2 ${gene_bed})
  end=\$(cut -f3 ${gene_bed})
  echo "${meta.gene} - chr: \$chr, start: \$start, end: \$end"

  # calculate coverage per base
  samtools depth \\
    -a --min-BQ ${params.min_BQ} --min-MQ ${params.min_MQ} \\
    -r \$chr:\$start-\$end ${bam} \\
  > ${meta.id}_${meta.gene}_coverage.tsv
  """
}

// workflow
workflow bulk {
  take:
  ch_bam_x_gene

  main:
  get_coverage(ch_bam_x_gene)

  emit:
  get_coverage.out
}