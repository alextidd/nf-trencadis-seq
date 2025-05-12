#!/bin/bash

# modules
module load samtools-1.19/python-3.12.0

# subset to test barcodes
bin/subset-bam \
  -b /lustre/scratch126/casm/team294rr/lm26/panel/mapped/TRAC-2-8179-Cell1/TRAC_2_8179_scisoseq.mapped.bam \
  --cell-barcodes test/data/cell_barcodes.txt \
  --out-bam test/data/PB_panel_KX004_subset.bam

# subset to genes
samtools view -L test/data/genes.bed -b test/data/PB_panel_KX004_subset.bam \
> test/data/PB_panel_KX004_subset.tmp.bam
mv test/data/PB_panel_KX004_subset.tmp.bam test/data/PB_panel_KX004_subset.bam

# index
samtools index test/data/PB_panel_KX004_subset.bam