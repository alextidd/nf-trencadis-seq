#!/bin/bash
# cd /nfs/casm/team268im/at31/nextflow/nf-trencadis-seq/test ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 02_run_tsunpo -o log/%J_02_run_tsunpo.out -e log/%J_02_run_tsunpo.err 'bash src/02_run_tsunpo.sh'

# dirs
tp_dir=/lustre/scratch125/casm/staging/team294/ty2/SSC/ngs/pacbio/nf-trencadis-seq/
work_dir=$LUSTRE_TEAM/nextflow/nf-trencadis-seq/test/work/

# run
nextflow run ../ \
  --samplesheet $tp_dir/samplesheet_trencadis.csv \
  --genes $tp_dir/genes_mn7_SMAD4.txt \
  --out_dir out/tsunpo/ \
  --location local \
  --refcds $tp_dir/data/refcds/GRCh38/RefCDS_human_GRCh38_GencodeV18_recommended.rda \
  --gencode_gff3 $tp_dir/data/gencode_gff3/GRCh38/gencode.v39.annotation.gff3.gz \
  -w $work_dir \
  -resume \
  -N at31@sanger.ac.uk