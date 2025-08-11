#!/bin/bash
# cd /nfs/casm/team268im/at31/nextflow/nf-trencadis-seq/test ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 03_run_tsunpo_large -o log/%J_03_run_tsunpo_large.out -e log/%J_03_run_tsunpo_large.err 'bash src/03_run_tsunpo_large.sh'

# dirs
tp_dir=/lustre/scratch125/casm/staging/team294/ty2/SSC/ngs/pacbio/nf-trencadis-seq/
work_dir=$LUSTRE_TEAM/nextflow/nf-trencadis-seq/test/work/

# run
nextflow run ../ \
  --samplesheet $tp_dir/samplesheet_trencadis_PD53623b.csv \
  --genes $tp_dir/genes_mn7.txt \
  --out_dir out/tsunpo_large/ \
  --location local \
  --refcds $tp_dir/data/refcds/GRCh38/RefCDS_human_GRCh38_GencodeV18_recommended.rda \
  --gencode_gff3 $tp_dir/data/gencode_gff3/GRCh38/gencode.v39.annotation.gff3.gz \
  -w $work_dir \
  -resume \
  -N at31@sanger.ac.uk