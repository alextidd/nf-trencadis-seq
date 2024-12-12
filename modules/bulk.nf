
// get coverage per site
process get_coverage {
  tag "${meta.id}_${meta.gene}"
  label 'week10gb'
  publishDir "${params.out_dir}/runs/${meta.id}/${meta.gene}/",
    mode: "copy",
    pattern: "*_coverage.tsv"

  input:
  tuple val(meta), path(bam), path(bai),
        path(mutations),
        path(celltypes), path(cell_barcodes),
        path(gene_bed), path(gene_features), path(exonic_bed)

  output:
  tuple val(meta), path("${meta.id}_${meta.gene}_coverage.tsv")

  script:
  """
  module load samtools-1.19/python-3.12.0

  # get gene info
  chr=\$(cut -f1 ${gene_bed})
  start=\$(cut -f2 ${gene_bed})
  end=\$(cut -f3 ${gene_bed})
  echo "${meta.gene} - chr: \$chr, start: \$start, end: \$end"

  # calculate coverage per base
  echo -e "chr\tpos\tcov" > ${meta.id}_${meta.gene}_coverage.tsv
  samtools depth \\
    -a --min-BQ ${params.min_BQ} --min-MQ ${params.min_MQ} \\
    -r \$chr:\$start-\$end ${bam} \\
  >> ${meta.id}_${meta.gene}_coverage.tsv
  """
}

// genotype mutations
process genotype_mutations {
  tag "${meta.id}_${meta.gene}"
  label 'normal10gb'
  publishDir "${params.out_dir}/runs/${meta.id}/${meta.gene}/",
    mode: "copy"
  
  input:
  tuple val(meta), path(bam), path(bai),
        path(mutations),
        path(celltypes), path(cell_barcodes),
        path(gene_bed), path(gene_features), path(exonic_bed)
  
  output:
  tuple val(meta), path("genotyped_mutations.tsv")
        
  script:
  """
  # genotype
  bulk_genotype_mutations.R \\
    --bam ${bam} \\
    --gene_bed ${gene_bed} \\
    --mutations ${mutations} \\
    --min_MQ ${params.min_MQ} \\
    --min_BQ ${params.min_BQ}
  """
}

// plot coverage
process plot_coverage {
  tag "${meta.gene}"
  label 'normal10gb'
  maxRetries 10
  publishDir "${params.out_dir}/genes/${meta.gene}/", mode: "copy"

  input:
  tuple val(meta), val(ids), path(cov)
  
  output:
  tuple val(meta), path("${meta.gene}_coverage.${params.plot_device}")
  
  script:
  """
  #!/usr/bin/env Rscript
  
  # libraries
  library(magrittr)
  library(ggplot2)

  # get coverage data
  files <- c("${cov.join('", "')}")
  ids <- c("${ids.join('", "')}")
  cov <-
    purrr::map2(files, ids, function(file, id) {
      readr::read_tsv(file) %>% dplyr::mutate(id = id)
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(lab = paste0(id, "\\nmedian ", median(cov),
                               " [", paste(range(cov), collapse = "-"), "]"))
                   

  # plot coverage
  p <-
    cov %>%
    ggplot(aes(x = pos, y = cov)) +
    geom_area() +
    theme_classic() +
    facet_grid(lab ~ .) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(strip.text.y.right = element_text(angle = 0),
          panel.spacing = unit(0, "lines"),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()) +
    ggtitle("${meta.gene} coverage")

  ggsave("${meta.gene}_coverage.${params.plot_device}", p,
    width = 10, height = 0.5 * dplyr::n_distinct(cov\$id), units = "in", dpi = 300)
  """
}

// workflow
workflow bulk {
  take:
  ch_bam_x_gene

  main:
  get_coverage(ch_bam_x_gene)
  genotype_mutations(ch_bam_x_gene)

  // plot coverage and mutations
  // ch_bam_x_gene.join(get_coverage.out).join(genotype_mutations.out).view()

  // plot coverage per gene
  get_coverage.out
  | map { meta, cov -> [meta.subMap('gene'), meta.id, cov] }
  | groupTuple()
  | plot_coverage
  
  emit:
  get_coverage.out
}