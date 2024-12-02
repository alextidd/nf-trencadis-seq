// subset the BAM to only reads from barcodes of interest
process subset_bam_to_barcodes {
  tag "${meta.id}_${meta.gene}"
  label 'normal'
  publishDir "${params.out_dir}/runs/${meta.id}/${meta.gene}/",
    mode: "copy",
    pattern: "*_subset.{bam,bam.bai}"

  input:
  tuple val(meta),
        path(bam), path(bai),
        path(celltypes),
        path(mutations),
        path(cell_barcodes),
        path(gene_bed), path(gene_regions), path(gene_features)

  output:
  tuple val(meta),
        path("${meta.id}_subset.bam"), path("${meta.id}_subset.bam.bai"),
        path(celltypes),
        path(mutations),
        path(cell_barcodes),
        path(gene_bed), path(gene_regions), path(gene_features)

  script:
  """
  module load samtools-1.19/python-3.12.0 

  # subset to cell barcodes
  subset-bam \
    -b ${bam} \
    --cell-barcodes ${cell_barcodes} \
    --out-bam ${meta.id}_subset.bam
  
  # index
  samtools index -@ ${task.cpus} ${meta.id}_subset.bam
  """
}

// split subsetted BAMs by cell, count coverage
process get_coverage_per_cell {
  tag "${meta.id}_${meta.gene}"
  label 'week10gb'
  publishDir "${params.out_dir}/runs/${meta.id}/${meta.gene}/",
    mode: "copy",
    pattern: "*_coverage_per_cell.tsv"

  input:
  tuple val(meta), path(bam), path(bai),
        path(celltypes),
        path(mutations),
        path(cell_barcodes),
        path(gene_bed), path(gene_regions), path(gene_features)

  output:
  tuple val(meta), path("${meta.id}_${meta.gene}_coverage_per_cell.tsv"),
        path(celltypes),
        path(cell_barcodes),
        path(gene_bed), path(gene_regions), path(gene_features),
        emit: cov
  tuple val(meta),
        path(celltypes), path(mutations), path(gene_bed),
        path("cell_bams/*"),
        emit: cell_bams

  script:
  """
  module load samtools-1.19/python-3.12.0

  # get gene info
  chr=\$(cut -f1 ${gene_bed})
  start=\$(cut -f2 ${gene_bed})
  end=\$(cut -f3 ${gene_bed})
  echo "${meta.gene} - chr: \$chr, start: \$start, end: \$end"

  echo "getting cell barcodes"
  # get unique cell barcodes in the BAM, trim CB:Z: prefix
  samtools view ${bam} | cut -f12- | tr "\\t" "\\n" |
  grep "CB:Z:" | sed 's/^CB:Z://g' | awk '!x[\$0]++' \
  > cell_barcodes.txt

  echo "creating a bam per cell"
  rm -rf cell_bams ; mkdir -p cell_bams
  while read -r CB ; do
    # subset
    echo \$CB > \$CB.txt
    subset-bam \
      -b ${bam} \
      --cell-barcodes \$CB.txt \
      --out-bam cell_bams/\$CB.bam
    # index
    samtools index -@ ${task.cpus} cell_bams/\$CB.bam
    rm \$CB.txt
  done < cell_barcodes.txt

  # initialise gene coverage file
  rm -rf ${meta.gene} ; mkdir ${meta.gene}
  touch ${meta.gene}_cov.tsv

  echo "calculating coverage per base per cell"
  while read -r CB ; do
    (
      echo "\$CB" ;
      samtools depth \\
        -a --min-BQ ${params.min_BQ} --min-MQ ${params.min_MQ} \\
        -r \$chr:\$start-\$end cell_bams/\$CB.bam \\
      | cut -f3 ;
    ) > ${meta.gene}/\$CB.tsv
  done < cell_barcodes.txt

  echo "creating coverage file"
  (
    echo -e "chr\\tpos\\tgene" ;
    echo -e "\$chr\\t\$start\\t\$end" \
    | awk -F"\\t" -v gene=${meta.gene} -v OFS="\\t" \
      '{for(i=\$2; i<=\$3; i++) print \$1, i, gene}' ;
  ) > ${meta.gene}_coords.tsv

  echo "adding coverage per cell to coverage file"
  # (sort for consistent ordering)
  find ${meta.gene} -name "*.tsv" | sort | xargs paste \\
  > ${meta.gene}_cov.tsv

  echo "combining outputs"
  paste ${meta.gene}_coords.tsv ${meta.gene}_cov.tsv > ${meta.gene}_cov_per_cell.tsv
  ls *_cov_per_cell.tsv | head -1 | xargs head -1 \\
  > ${meta.id}_${meta.gene}_coverage_per_cell.tsv
  cat *_cov_per_cell.tsv | grep -vP "^chr\\tpos" \\
  >> ${meta.id}_${meta.gene}_coverage_per_cell.tsv
  """
}

// summarise coverage per celltype
process get_coverage_per_celltype {
  tag "${meta.id}_${meta.gene}"
  label 'normal20gb'
  maxRetries 10
  
  input:
  tuple val(meta),
        path(cov),
        path(celltypes),
        path(cell_barcodes),
        path(gene_bed), path(gene_regions), path(gene_features)

  output:
  tuple val(meta),
        path(gene_bed), path(gene_regions), path(gene_features),
        path("${meta.id}_${meta.gene}_cov_per_ct.rds"),
        emit: cov_per_ct
  tuple val(meta),
        path("${meta.id}_${meta.gene}_coverage_per_exonic_position.tsv"),
        emit: exonic_coverage

  script:
  """
  get_coverage_per_celltype.R \\
    --gene ${meta.gene} \\
    --regions ${gene_regions} \\
    --cov ${cov} \\
    --cell_barcodes ${cell_barcodes} \\
    --celltypes ${celltypes} \\
    --meta_id ${meta.id}
  """
}

// workflow
workflow single_cell {
  take:
  ch_bam_x_gene

  main:
  // get per-cell bams
  subset_bam_to_barcodes(ch_bam_x_gene)

  // get coverage per cell
  get_coverage_per_cell(subset_bam_to_barcodes.out)

  // summarise coverage
  get_coverage_per_celltype(get_coverage_per_cell.out.cov)
  
  emit:
  cov_per_ct = get_coverage_per_celltype.out.cov_per_ct
  cell_bams = get_coverage_per_cell.out.cell_bams
}