#!/usr/bin/env nextflow
// TODO: add genotyping step to pipeline
// TODO: simplify mutations + fix handling in plot_mut_pie()
// TODO: ignore noncoding mutations!
// TODO: fold Rmd reports into the nextflow pipeline
// TODO: fix called mutation annotations to fit + facet by celltype
// TODO: add gene expression distribution alongside mutations / coverage plots
// TODO: create summary plot per gene
// TODO: add celltype marker genes to the expression dotplot (ask Laura)
// TODO: check over the transcript IDs that Carol from pipeline sent - begin deconvolution in the dataset

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// Download a given sample's BAM from iRODS
// Then either retrieve the BAI or make one via indexing
// The maxForks of 10 was set after asking jc18 about best iRODS practices
process irods {
  tag "${meta.run}"
  maxForks 10
  label 'normal4core'

  input:
  tuple val(meta), val(bam),
        path(celltypes), path(mutations), path(cell_barcodes)
  
  output:
  tuple val(meta), path("${meta.run}.bam"), path("${meta.run}.bam.bai"),
        path(celltypes), path(mutations), path(cell_barcodes)
  
  script:
  """
  iget -K ${bam} ${meta.run}.bam
  if [[ `ils ${bam}.bai | wc -l` == 1 ]]
  then
      iget -K ${bam}.bai ${meta.run}.bam.bai
  else
      samtools index -@ ${task.cpus} ${meta.run}.bam
  fi
  """
}

// The equivalent of an irods download, but for a local copy of samplesheet
// Symlink the BAM/BAI appropriately so they're named the right thing for downstream
process local {
  tag "${meta.run}"
  maxForks 10
  label 'normal4core'

  input:
  tuple val(meta), val(bam),
        path(celltypes), path(mutations), path(cell_barcodes)

  output:
  tuple val(meta), path("${meta.run}.bam"), path("${meta.run}.bam.bai"),
        path(celltypes), path(mutations), path(cell_barcodes)

  script:
  """
  # create local symbolic link 
  ln -s ${bam} ${meta.run}.bam
  if [ -f "${bam}.bai" ] ; then
      ln -s ${bam}.bai ${meta.run}.bam.bai
  else
      samtools index -@ ${task.cpus} ${meta.run}.bam
  fi
  """
}

// check bam is not truncated before proceeding
process check_bam {
  tag "${meta.run}"
  label 'normal'
  errorStrategy 'ignore'

  input:
  tuple val(meta), path(bam), path(bai),
        path(celltypes),
        path(mutations),
        path(cell_barcodes)

  output:
  tuple val(meta), path(bam), path(bai),
        path(celltypes),
        path(mutations),
        path(cell_barcodes)

  script:
  """
  module load samtools-1.19/python-3.12.0 
  samtools quickcheck ${bam}
  """
}

// get gene coords
process get_gene_coords {
  tag "${gene}"
  label 'normal10gb'
  publishDir "${params.out_dir}/genes/${gene}/",
    mode: "copy",
    pattern: "*.bed"

  input:
  val(gene)
  path(gencode_gff3)
  path(refcds)

  output:
  tuple val(gene),
        path("${gene}.bed"),
        path("${gene}_features.tsv"),
        path("${gene}_exonic_regions.bed"), emit: gene_features
  path "${gene}.bed", emit: gene_bed

  script:
  """
  #!/usr/bin/env Rscript

  # libraries
  library(GenomicRanges)
  library(magrittr)

  # load the gff
  gff <- ape::read.gff("${gencode_gff3}")

  # load the refcds
  load("${refcds}")
  g <- RefCDS[[which(purrr::map_lgl(RefCDS, ~ .x\$gene_name == "${gene}"))]]
  ccds <-
    tibble::tibble(chr = paste0("chr", g\$chr), start = g\$intervals_cds[, 1],
                   end = g\$intervals_cds[, 2], type = "ccds")

  # get feature coords
  features <-
    gff %>%
    dplyr::filter(
      grepl("gene_name=${gene};", attributes),
      # get protein-coding gene coords
      (type == "gene" & grepl("gene_type=protein_coding;", attributes)) |
      # get transcript / exon coords
      (type %in% c("transcript", "exon"))) %>%
    dplyr::transmute(chr = seqid, start, end, gene = "${gene}", strand,
                     type, attributes) %>%
    dplyr::bind_rows(ccds)
  
  # get the gene coords
  gene_bed <-
    features %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(chr, start, end, gene, strand)
  
  # collapse all exonic regions
  exonic_regions <- 
    features %>%
    dplyr::filter(type == "exon") %>%
    {GRanges(seqnames = .\$chr,
             ranges = IRanges(.\$start, .\$end))} %>%
    reduce() %>%
    {tibble::tibble(chr = as.vector(seqnames(.)),
                    start = start(.),
                    end = end(.))}

  # write
  readr::write_tsv(features, "${gene}_features.tsv")
  readr::write_tsv(gene_bed, "${gene}.bed", col_names = FALSE)
  readr::write_tsv(exonic_regions, "${gene}_exonic_regions.bed", col_names = FALSE)
  """
}

// collate gene coords
process collate_gene_coords {
  label 'normal'
  publishDir "${params.out_dir}/genes/", mode: "copy"

  input:
  path(genes)

  output:
  path("genes.bed")

  script:
  """
  echo -e "chr\tstart\tend\tgene\tstrand" > genes.bed
  cat ${genes} >> genes.bed
  """
}

// subset the BAM to only reads in the genes
process subset_bam_to_genes {
  tag "${meta.run}_${meta.gene}"
  label 'normal'

  input:
  tuple val(meta), path(bam), path(bai),
        path(celltypes),
        path(mutations),
        path(cell_barcodes),
        path(gene_bed), path(gene_features), path(exonic_bed)

  output:
  tuple val(meta),
        path("${meta.run}_subset_genes.bam"),
        path("${meta.run}_subset_genes.bam.bai"),
        path(celltypes),
        path(mutations),
        path(cell_barcodes),
        path(gene_bed), path(gene_features), path(exonic_bed)

  script:
  """
  module load samtools-1.19/python-3.12.0 

  # subset to reads that fall in gene regions
  samtools view -L ${gene_bed} -b ${bam} > ${meta.run}_subset_genes.bam

  # index
  samtools index -@ ${task.cpus} ${meta.run}_subset_genes.bam
  """
}

// subset the BAM to only reads from barcodes of interest
process subset_bam_to_barcodes {
  tag "${meta.run}_${meta.gene}"
  label 'normal'
  publishDir "${params.out_dir}/runs/${meta.run}/${meta.gene}/",
    mode: "copy",
    pattern: "*_subset.{bam,bam.bai}"

  input:
  tuple val(meta),
        path(bam), path(bai),
        path(celltypes),
        path(mutations),
        path(cell_barcodes),
        path(gene_bed), path(gene_features), path(exonic_bed)

  output:
  tuple val(meta),
        path("${meta.run}_subset.bam"), path("${meta.run}_subset.bam.bai"),
        path(celltypes),
        path(mutations),
        path(cell_barcodes),
        path(gene_bed), path(gene_features), path(exonic_bed)

  script:
  """
  module load samtools-1.19/python-3.12.0 

  # subset to cell barcodes
  subset-bam \
    -b ${bam} \
    --cell-barcodes ${cell_barcodes} \
    --out-bam ${meta.run}_subset.bam
  
  # index
  samtools index -@ ${task.cpus} ${meta.run}_subset.bam
  """
}

// split subsetted BAMs by cell, count coverage
process get_coverage_per_cell {
  tag "${meta.run}_${meta.gene}"
  label 'week10gb'
  publishDir "${params.out_dir}/runs/${meta.run}/${meta.gene}/",
    mode: "copy",
    pattern: "*_coverage_per_cell.tsv"

  input:
  tuple val(meta), path(bam), path(bai),
        path(celltypes),
        path(mutations),
        path(cell_barcodes),
        path(gene_bed), path(gene_features), path(exonic_bed)

  output:
  tuple val(meta), path("${meta.run}_${meta.gene}_coverage_per_cell.tsv"),
        path(celltypes),
        path(cell_barcodes),
        path(gene_bed), path(gene_features), path(exonic_bed),
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
  > ${meta.run}_${meta.gene}_coverage_per_cell.tsv
  cat *_cov_per_cell.tsv | grep -vP "^chr\\tpos" \\
  >> ${meta.run}_${meta.gene}_coverage_per_cell.tsv
  """
}

// genotype the mutations
process genotype_mutations {
  tag "${meta.run}_${meta.gene}"
  label 'normal10gb'
  publishDir "${params.out_dir}/runs/${meta.run}/${meta.gene}/",
    mode: "copy"
  
  input:
  tuple val(meta),
        path(celltypes), path(mutations), path(gene_bed),
        path(cell_bams, stageAs: "cell_bams/*")
  
  output:
  tuple val(meta),
        path("genotyped_mutations_per_celltype.tsv"), emit: geno_per_ct
  path("genotyped_mutations_per_cell.tsv"), emit: geno_per_cell
  path("mutations.tsv"), emit: muts
        
  
  script:
  """
  # genotype
  genotype_mutations.R \\
    --gene_bed ${gene_bed} \\
    --mutations ${mutations} \\
    --celltypes ${celltypes} \\
    --min_MQ ${params.min_MQ} \\
    --min_BQ ${params.min_BQ}
  """
}

// wrangle data for plots
process wrangle_data {
  tag "${meta.run}_${meta.gene}"
  label 'normal20gb'
  maxRetries 10
  
  input:
  tuple val(meta),
        path(cov),
        path(celltypes),
        path(cell_barcodes),
        path(gene_bed), path(gene_features), path(exonic_bed),
        path(geno_per_ct)
  path(refcds)

  output:
  tuple val(meta), path("*.rds"), path(geno_per_ct),
        emit: plot_data
  tuple val(meta),
        path("${meta.run}_${meta.gene}_coverage_per_exonic_position.tsv"),
        emit: exonic_coverage

  script:
  """
  wrangle_data.R \\
    --gene ${meta.gene} \\
    --gene_bed ${gene_bed} \\
    --exonic_bed ${exonic_bed} \\
    --refcds ${refcds} \\
    --gene_features ${gene_features} \\
    --cov ${cov} \\
    --cell_barcodes ${cell_barcodes} \\
    --celltypes ${celltypes} \\
    --meta_run ${meta.run}
  """
}

// plot coverage
process plot_coverage {
  tag "${meta.run}_${meta.gene}"
  label 'normal10gb'
  maxRetries 10
  publishDir "${params.out_dir}/runs/${meta.run}/${meta.gene}/", mode: "copy"
  errorStrategy 'ignore'

  input:
  tuple val(meta), path(rdss), path(geno_per_ct)
  
  output:
  path("${meta.run}_${meta.gene}_*_plot.${params.plot_device}")
  path("${meta.run}_${meta.gene}_pie_mutations_plot.rds"), optional: true

  script:
  """
  plot_coverage.R \\
    --meta_run ${meta.run} \\
    --gene ${meta.gene} \\
    --geno_per_ct ${geno_per_ct} \\
    --min_cov ${params.min_cov} \\
    --plot_device ${params.plot_device}
  """
}

process plot_5_prime_dropoff {
  tag "${gene}"
  label 'normal10gb'
  maxRetries 10
  publishDir "${params.out_dir}/genes/${gene}/",
    mode: "copy",
    pattern: "*.{png,rds}"
  errorStrategy 'retry'

  input:
  tuple val(gene),
        path(exonic_coverages)

  output:
  path("${gene}_exonic_coverage_plot.${params.plot_device}")
  path("${gene}_exonic_coverage_plot.rds")

  script:
  """
  plot_5_prime_dropoff.R \\
    --gene ${gene} \\
    --min_cov ${params.min_cov} \\
    --plot_device ${params.plot_device}
  """
}

// main workflow
workflow {
  
  // validate input parameters
  validateParameters()

  // print summary of supplied parameters
  log.info paramsSummaryLog(workflow)

  // get metadata + bam paths  
  Channel.fromPath(params.samplesheet, checkIfExists: true)
  | splitCsv(header: true)
  | map { row ->
      def run = row.kit + "_" + row.seq_type + "_" + row.id
      def meta = [id:row.id, seq_type:row.seq_type, kit:row.kit, run:run]
      [meta,
       file(row.bam, checkIfExists: true),
       file(row.celltypes, checkIfExists: true),
       file(row.mutations, checkIfExists: true),
       file(row.cell_barcodes, checkIfExists: true)]
  }
  | set { input }
  
  // download or locally link bams
  if (params.location == "irods") {
    // download input from irods
    ch_bam = irods(input)
  }
  else if (params.location == "local") {
    // input are locally available
    ch_bam = local(input)
  }

  // check bam is not truncated before proceeding
  ch_checked_bam = check_bam(ch_bam)
  
  // get gencode gff3
  gencode_gff3 = file(params.gencode_gff3)

  // get refcds
  refcds = file(params.refcds)

  // get gene coords
  ch_gene = Channel.fromPath(params.genes).splitText().map{it -> it.trim()}
  get_gene_coords(ch_gene, gencode_gff3, refcds)

  // output bed of all gene coords
  get_gene_coords.out.gene_bed.collect()
  | collate_gene_coords

  // combine gene-x-id, subset bams
  ch_checked_bam.combine(get_gene_coords.out.gene_features)
  | map { meta, bam, bai, celltypes, mutations, cell_barcodes,
          gene, gene_bed, gene_features, exonic_bed ->
          [meta + [gene: gene],
           bam, bai, celltypes, mutations, cell_barcodes, 
           gene_bed, gene_features, exonic_bed]
  }
  | subset_bam_to_genes
  | subset_bam_to_barcodes

  // get coverage per cell
  get_coverage_per_cell(subset_bam_to_barcodes.out)

  // genotype cells
  genotype_mutations(get_coverage_per_cell.out.cell_bams)
  
  // prep data for plotting
  wrangle_data(
    get_coverage_per_cell.out.cov.join(genotype_mutations.out.geno_per_ct),
    refcds)

  // plot coverage
  plot_coverage(wrangle_data.out.plot_data)

  // // plot coverage and mutations
  // wrangle_data(get_coverage_per_cell.out.cov)
  // plot_coverage(wrangle_data.out.plot_data)

  // // plot 5' drop-off per gene
  // ch_exonic_cov_per_gene = wrangle_data.out.exonic_coverage.groupTuple()
  // plot_5_prime_dropoff(ch_exonic_cov_per_gene)
  
}