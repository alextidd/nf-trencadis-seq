# nf-trencadis-seq

## Overview

`nf-trencadis-seq` is a Nextflow pipeline designed for performing genotyping and 
coverage analyses for the Trencadis-Seq protocol, which combines parallel 10X
short-read and PacBio long-read sequencing of the same barcoded library to 
genotype mutations in single cells with associated celltype information.

1. If `--location irods`, retrieve the BAM from iRODS.
2. Index the BAM.
3. Get the coordinates of the genes of interest.
4. Subset the BAM to only reads overlapping the genes of interest.
5. Subset the BAM to only the barcodes of interest.
6. Get coverage per cell and per celltype at every site in the genes of interest.
7. Genotype all mutations per cell and per celltype in the genes of interest.
8. Plot coverage and genotyping results for each gene in each sample.
9. Knit the report.

## Usage

### Samplesheet

First, prepare a comma-delimited samplesheet with your input data. It should 
look like this:

| id             | bam                                                                                                              | cell_barcodes                                                                                                     | celltypes                                                                                                  | mutations                                                                                           |
|---------------|----------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------|
| TX_WT_AX001  | /lustre/scratch126/casm/team294rr/lm26/Laura-10X_genomics_scRNA/cellranger720_count_48211_7323STDY14579485_GRCh38-2020-A/possorted_genome_bam.bam  | /lustre/scratch125/casm/team268im/at31/trencadis-seq/out/blood/celltype_annotations/TX_WT_AX001_cell_barcodes_seurat.txt  | /lustre/scratch125/casm/team268im/at31/trencadis-seq/out/blood/celltype_annotations/AX001_barcode_TX_celltype_annotations.csv  | /lustre/scratch125/casm/team268im/at31/trencadis-seq/out/blood/mutations/AX001_mutations.tsv  |
| PB_WT_AX001  | /lustre/scratch126/casm/team294rr/lm26/whole_transcriptome/mapped/TRAC-2-8004-Cell1/TRAC_2_8004_scisoseq.mapped.bam  | /lustre/scratch125/casm/team268im/at31/trencadis-seq/out/blood/celltype_annotations/PB_WT_AX001_cell_barcodes_seurat.txt  | /lustre/scratch125/casm/team268im/at31/trencadis-seq/out/blood/celltype_annotations/AX001_barcode_PB_celltype_annotations.csv  | /lustre/scratch125/casm/team268im/at31/trencadis-seq/out/blood/mutations/AX001_mutations.tsv  |
| PB_panel_AX001  | /lustre/scratch126/casm/team294rr/lm26/panel/mapped/TRAC-2-8178-Cell1/TRAC_2_8178_scisoseq.mapped.bam  | /lustre/scratch125/casm/team268im/at31/trencadis-seq/out/blood/celltype_annotations/PB_panel_AX001_cell_barcodes_seurat.txt  | /lustre/scratch125/casm/team268im/at31/trencadis-seq/out/blood/celltype_annotations/AX001_barcode_PB_celltype_annotations.csv  | /lustre/scratch125/casm/team268im/at31/trencadis-seq/out/blood/mutations/AX001_mutations.tsv  |

Each row represents a BAM file corresponding to a unique ID, with associated 
celltype annotations, barcodes, and mutations to genotype. The samplesheet must
contain at least the columns `id`, `bam`, and `mutations`.

#### Mutations

The mutations files should be tab-delimited and look like this:

```
chr     pos     ref     alt
chr2    43225731        G       A
chr20   32369041        A       G
chr20   32433398        C       T
chr20   32433398        C       T
chr20   32433398        C       T
chr22   28694073        G       A
chr22   28710054        TG      T
```

#### Celltypes

The celltypes files should be comma-delimited and look like this:

```
barcode,celltype
AGAGAGGCTTCCACAA,NK_cell
GTCGGGTACAAATGGT,NK_cell
TGCTCTCGAAATACCC,NK_cell
GTTATTGCTCATCGGC,T_cells
CTTAATCTGGGACGTA,T_cells
AAACTGCTGAGATTGA,B_cell
```

Make sure the barcodes in this file match those present in the cellltypes file
and in the BAM (check the `CB:Z` tag of the reads in the BAM to compare). 10X
BAMs will have a `-1` suffix on the barcodes, so if you are analysing 10X BAMs,
be sure to include this.

**N.B.** If you do not have celltype annotations available, you can provide a
samplesheet without the `celltypes` column and pass the parameter
`--no_celltypes`. The pipeline will group all cells together and treat them as
one celltype, giving coverage and genotyping information for each gene across
the entire population.

#### Cell barcodes

The cell barcodes file should be a text file listing all barcodes of interest:

```
AGCTACCCTTCGTAGT
TTTGGGAACCATCGAT
GATGCAACTCGTGACC
ACTCATGGAAAATTGG
ACATAGAGAGAAGTAT
```

Before genotyping and coverage analyses, the pipeline will subset the BAM to 
only those barcodes of interest. You can use this if you would like to keep
only high-quality cells that passed QC or only cells with annotations, for
example.

**N.B.** If you do not have cell barcodes of interest, you can provide a 
samplesheet without the `cell_barcodes` column and pass the parameter
`--no_cell_barcodes`, and the pipeline will skip this subsetting step.

### Genes

The other required input file is the genes of interest. These are the genes in
which you would like to assess coverage and genotype mutations across samples.
The genes file must be a text file with one gene per line, like this:

```
DNMT3A
TET2
SF3B1
SRSF2
TP53
```

### Run

Now, you can run the pipeline using:

```bash
nextflow run ${wd}/../nextflow/nf-trencadis-seq \
  --samplesheet samplesheet.csv \
  --genes driver_genes.txt
```

To get a full list of parameters, use the `--help` flag:

```bash
$ nextflow run . --help
```

## Parameters

### Required options

- `--samplesheet` [string]: Comma-separated file containing the columns 'id', 'bam', 'celltypes', 'mutations', and 'cell_barcodes'.
- `--genes` [string]: Text file of genes of interest, with one gene per line.

### Input/output options

- `--location` [string]: Are the BAMs saved locally or on iRODs? [default: local]
- `--out_dir` [string]: Output directory. [default: out/]
- `--plot_device` [string]: File type for saving plots. (accepted: pdf, png) [default: png]
- `--min_cov` [integer]: Minimum coverage of a site in a cell for it to be genotyped. [default: 1]
- `--min_MQ` [integer]: Minimum mapping quality. [default: 0]
- `--min_BQ` [integer]: Minimum base quality. [default: 0]
- `--no_celltypes` [boolean]: No celltype annotations in the samplesheet (`celltypes` column).
- `--no_cell_barcodes` [boolean]: No cell barcodes in the samplesheet (`cell_barcodes` column).

### Reference files

- `--no_chr` [boolean]: Is there a 'chr' prefix on the chromosome names? [default: false]
- `--refcds` [string]: Path to RefCDS Rda object from the dndscv package. [default: ${baseDir}/data/refcds/GRCh38/RefCDS_human_GRCh38_GencodeV18_recommended.rda] 
- `--gencode_gff3` [string]: Path to a GENCODE GFF3 annotation file. [default: ${baseDir}/data/gencode_gff3/GRCh38/gencode.v39.annotation.gff3.gz] 

If you wish to run the pipeline using a different genome build, you will have 
to change these. 

**N.B.** Make sure the contig names (e.g. chr1 vs 1) are consistent across the
BAM, the mutations file, the RefCDS, and the GENCODE GFF3.

### Generic options

- `--help` [boolean]: Display help text.

## Output

An example of an output directory, for the sample `PB_panel_AX001` and the gene
`DNMT3A`:

```
.
├── report.html
├── genes
│   └── DNMT3A
│       └── DNMT3A.bed
└── runs
    └── PB_panel_KX004
        └── DNMT3A
            ├── PB_panel_KX004_DNMT3A_subset.bam.bai
            ├── PB_panel_KX004_DNMT3A_subset.bam
            ├── PB_panel_KX004_DNMT3A_coverage_per_cell.tsv
            ├── PB_panel_KX004_DNMT3A_coverage_per_celltype.tsv
            ├── mutations.tsv
            ├── PB_panel_KX004_DNMT3A_genotyped_mutations_per_cell.tsv
            ├── PB_panel_KX004_DNMT3A_genotyped_mutations_per_celltype.tsv
            └── *_plot.png
```

- `report.html` - Summary of all plotting outputs.
- `DNMT3A.bed` - Gene coordinates used by the pipeline.
- `PB_panel_KX004_DNMT3A_subset.bam` - BAM, subset to reads overlapping the
gene and with barcodes present in the cell barcodes file.
- `PB_panel_KX004_DNMT3A_coverage_per_cell.tsv` - Coverage per cell per site.
- `PB_panel_KX004_DNMT3A_coverage_per_celltype.tsv` - Coverage per celltype per
site.
- `mutations.tsv` - All mutations in the gene that were genotyped.
- `PB_panel_KX004_DNMT3A_genotyped_mutations_per_cell.tsv` - Reference and 
alternative reads for each mutation in the gene per cell.
- `PB_panel_KX004_DNMT3A_genotyped_mutations_per_celltype.tsv` - Reference and 
alternative reads for each mutation in the gene per celltype.
- `*_plot.png` - Various plots of the genotyping and coverage information.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.