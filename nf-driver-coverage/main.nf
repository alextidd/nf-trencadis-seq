#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// validate input parameters
if (params.validate_params) {
    validateParameters()
}

// Download a given sample's BAM from iRODS
// Then either retrieve the BAI or make one via indexing
// The maxForks of 10 was set after asking jc18 about best iRODS practices
process irods {
  tag "${meta.id}"
  maxForks 10
  label 'normal4core'
  input:
    tuple val(meta), val(bam)
  output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bams
  script:
    """
    iget -K ${bam} ${meta.id}.bam
    if [[ `ils ${bam}.bai | wc -l` == 1 ]]
    then
        iget -K ${bam}.bai ${meta.id}.bam.bai
    else
        samtools index -@ ${task.cpus} ${meta.id}.bam
    fi
    """
}

// The equivalent of an irods download, but for a local copy of samplesheet
// Symlink the BAM/BAI appropriately so they're named the right thing for downstream
process local {
  tag "${meta.id}"
  maxForks 10
  label 'normal4core'
  errorStrategy = {task.attempt <= 1 ? 'retry' : 'ignore'}
  input:
    tuple val(meta), val(bam)
  output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai")
  script:
    """
    # create local symbolic link 
    ln -s ${bam} ${meta.id}.bam
    if [ -f "${bam}.bai" ] ; then
        ln -s ${bam}.bai ${meta.id}.bam.bai
    else
        samtools index -@ ${task.cpus} ${meta.id}.bam
    fi
    """
}

// get driver gene coords
process get_driver_gene_coords {
  label 'normal'
  input:
    path(drivers_txt)
  output:
    path("driver_genes.bed")
  script:
    """
    (
      while read -r gene ; do

      # get the gene coords from the gff3 file
      # (chr start end gene strand)
      zcat /lustre/scratch125/casm/team268im/at31/reference/gencode/gencode.v28.annotation.gff3.gz \
      | grep -w \$gene \
      | awk -F"\t" -v OFS="\\t" -v gene=\$gene \
        '\$3 == "gene" {print \$1, \$4, \$5, gene, \$7}' ;

      done < $drivers_txt
    ) | cat > driver_genes.bed
    """
}

// subset the BAMs to only the driver regions
process subset_bams_to_drivers {
  tag "${meta.id}"
  label 'normal'
  input:
    tuple val(meta), path(bam), path(bai)
    path(drivers_bed)
  output:
    tuple val(meta),
          path("${meta.id}_subset.bam"), path("${meta.id}_subset.bam.bai")
    path(drivers_bed)
  script:
    """
    module load samtools-1.19/python-3.12.0 
    samtools view -L $drivers_bed -b $bam > ${meta.id}_subset.bam
    samtools index -@ ${task.cpus} ${meta.id}_subset.bam
    """
}

// split subsetted BAMs by cell, count coverage
process get_coverage_per_cell {
  tag "${meta.id}"
  label 'normal'
  publishDir "${params.out_dir}/", mode: "copy"
  errorStrategy 'ignore'
  input:
    tuple val(meta), path(bam), path(bai)
    path(drivers_bed)
  output:
    tuple val(meta), path("${meta.id}_driver_coverage_per_cell.tsv")
    path(drivers_bed)
  script:
    """
    module load samtools-1.19/python-3.12.0 

    echo "getting cell barcodes"
    mkdir cell_bams
    samtools view $bam | cut -f12 | cut -d':' -f3 \
    | awk '!x[\$0]++' \
    > cell_barcodes.txt

    echo "creating a bam per cell"
    while read -r CB ; do
      # subset
      echo \$CB > \$CB.txt
      subset-bam \
        -b $bam \
        --cell-barcodes \$CB.txt \
        --out-bam cell_bams/\$CB.bam
      # index
      samtools index -@ ${task.cpus} cell_bams/\$CB.bam
      rm \$CB.txt
    done < cell_barcodes.txt

    # run per gene
    while read -r chr start end gene strand ; do

      echo "processing gene \$gene"

      # initialise gene coverage file
      mkdir \${gene}
      touch \${gene}_cov.tsv

      echo "calculating coverage per base per cell"
      while read -r CB ; do
        (
          echo "\$CB" ;
          samtools depth -a -r \$chr:\$start-\$end cell_bams/\$CB.bam \
          | cut -f3 ;
        ) > \${gene}/\$CB.tsv
      done < cell_barcodes.txt

      echo "creating coverage file"
      (
        echo -e "chr\\tpos\\tgene" ;
        echo -e "\$chr\\t\$start\\t\$end" \
        | awk -F"\\t" -v gene=\$gene -v OFS="\\t" \
          '{for(i=\$2; i<=\$3; i++) print \$1, i, gene}' ;
      ) > \${gene}_coords.tsv

      echo "adding coverage per cell to coverage file"
      # (sort for consistent ordering)
      paste \$(ls \${gene}/*.tsv | sort) > \${gene}_cov.tsv

      # combine coords and coverage
      paste \${gene}_coords.tsv \${gene}_cov.tsv > \${gene}_cov_per_cell.tsv

    done < $drivers_bed

    echo "combining outputs"
    ls *_cov_per_cell.tsv | head -1 | xargs head -1 \
    > ${meta.id}_driver_coverage_per_cell.tsv
    cat *_cov_per_cell.tsv | grep -vP "^chr\\tpos" \
    >> ${meta.id}_driver_coverage_per_cell.tsv
    """
}

// plot coverage
process plot_coverage {
  tag "${meta.id}"
  label 'normal'
  publishDir "${params.out_dir}/", mode: "copy"
  beforeScript "module load R/4.2.2"
  input:
    tuple val(meta), path(cov)
    path(drivers_bed)
  output:
    path("*_coverage_plot.png")
  script:
    """
    #!/usr/bin/env Rscript

    # libraries
    library(magrittr)
    library(ggplot2)

    # get gene strand info (to orient 5'/3' ends)
    drivers <- 
      readr::read_tsv("${drivers_bed}",
                      col_names = c("chr", "start", "end", "gene", "strand"))

    # get coverage data
    cov <- readr::read_tsv("${cov}")

    # get celltype annotations
    annot <- readr::read_tsv("${params.celltypes}") %>%
      dplyr::filter(donor_id == "${meta.id}") %>%
      dplyr::select(CB = barcode_revc, celltype)

    # prep plotting data
    p_dat <-
      cov %>%
      tidyr::pivot_longer(cols = -c(chr, pos, gene),
                          names_to = "CB", values_to = "cov per cell") %>%
      dplyr::inner_join(annot)
    
    p_dat %>%
      {split(., .\$gene)} %>%
      purrr::map(function(p_dat_g) {
        gene <- unique(p_dat_g\$gene)

        # plot coverage per base
        p_cov_per_base <-
          p_dat_g %>%
          ggplot(aes(x = start, fill = `cov per cell`, group = `cov per cell`)) +
          geom_bar(width = 1) +
          viridis::scale_fill_viridis() +
          facet_grid(celltype ~ ., space = "free_y", scales = "free_y") +
          theme_minimal() +
          coord_cartesian(xlim = c(min(p_dat\$start), max(p_dat\$start)))

        # plot genotyping efficiency (% of cells with >= `min_cov` reads per position)
        genotyping_efficiency <-
          cov %>%
          dplyr::select(-c(chr, pos, gene)) %>%
          as.matrix() %>%
          {rowMeans(. >= min_cov)}
        p_genotyping_efficiency <-
          cbind(cov[, c("chr", "pos", "gene")], genotyping_efficiency) %>%
          ggplot(aes(x = pos, fill = genotyping_efficiency, y = 1)) +
          geom_col(width = 1) +
          theme_minimal() +
          coord_cartesian(xlim = c(min(p_dat\$pos), max(p_dat\$pos)))

        png(paste0("${meta.id}_", gene, "_coverage_plot.png"))
        cowplot::plot_grid(p_cov_per_base, p_genotyping_efficiency, ncol = 1,
                          rel_heights = c(4, 1), align = "v", axis = "b")
        dev.off()
      })
    """
}

// main workflow
workflow {
  
  // get metadata + bam paths  
  Channel.fromPath(params.samplesheet, checkIfExists: true)
  | splitCsv(header: true)
  | map { row ->
      meta = row.subMap('id') 
      [meta, file(row.bam, checkIfExists: true)]
  }
  | set { samplesheet }
  
  // download or locally link bams
  if (params.location == "irods") {
    // download samplesheet from irods
    bams = irods(samplesheet)
  }
  else if (params.location == "local") {
    // samplesheet are locally available
    sample_bams = local(samplesheet)
  }
  
  // get driver gene coords
  // TODO: create channel of genes, parallelise
  // Channel.fromPath(params.drivers).splitText().set { ch_drivers }
  // ch_drivers.view()

  drivers_txt = file(params.drivers)
  drivers_bed = get_driver_gene_coords(drivers_txt)

  // plot coverage per base per gene per cell
  subset_bams_to_drivers(sample_bams, drivers_bed)
  get_coverage_per_cell(subset_bams_to_drivers.out)
  plot_coverage(get_coverage_per_cell.out)
  
}
