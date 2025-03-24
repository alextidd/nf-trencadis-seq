# nf-trencadis-seq

## Overview

`nf-trencadis-seq` is a Nextflow pipeline designed for processing and analyzing sequencing data. This pipeline supports various input/output options, reference files, and generic options to customize the workflow according to your needs.

## Usage

To get a full list of parameters, use the `--help` flag:

```bash
$ nextflow run . --help
```

## Parameters

### Input/output options

- `--samplesheet` [string]: Comma-separated file containing the columns 'id', 'bam', 'celltypes', 'mutations', and 'cell_barcodes'.
- `--location` [string]: Are the BAMs saved locally or on iRODs? [default: local]
- `--out_dir` [string]: Output directory. [default: out/]
- `--genes` [string]: Text file of genes of interest, with one gene per line.
- `--plot_device` [string]: File type for saving plots. (accepted: pdf, png) [default: png]
- `--min_cov` [integer]: Minimum coverage of a site in a cell for it to be genotyped. [default: 1]
- `--min_MQ` [integer]: Minimum mapping quality. [default: 0]
- `--min_BQ` [integer]: Minimum base quality. [default: 0]
- `--no_celltypes` [boolean]: No celltype annotations in the samplesheet (`celltypes` column).
- `--no_cell_barcodes` [boolean]: No cell barcodes in the samplesheet (`cell_barcodes` column).
- `--no_chr` [boolean]: Is there a 'chr' prefix on the chromosome names?

### Reference files

- `--refcds` [string]: Path to RefCDS Rda object from the dndscv package.
- `--gencode_gff3` [string]: Path to a GENCODE GFF3 annotation file.

### Generic options

- `--help` [boolean]: Display help text.

## Directory Structure

The repository has the following structure:

```
main.nf
modules.json
nextflow_schema.json
nextflow.config
README.md
assets/
    report.Rmd
    schema_samplesheet.json
bin/
    bulk_genotype_mutations.R
    genotype_mutations.R
    get_coverage_per_celltype.R
    get_gene_coords.R
    plot_5_prime_dropoff.R
    plot_coverage_and_genotyping.R
    subset-bam
modules/
    local/
        ...
    nf-core/
        ...
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.