#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/nextflow/nf-trencadis-seq ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 01_run_test -o test/log/%J_01_run_test.out -e test/log/%J_01_run_test.err 'bash test/src/01_run_test.sh'

# dir
wd=$(pwd)
ref_dir=${wd}/../../reference

# run test
(
  cd test/out/
  nextflow run ../.. \
    --samplesheet samplesheet.csv \
    --out_dir ./ \
    --genes ../data/genes.txt \
    --gencode_gff3 $ref_dir/gencode/GRCh38/gencode.v39.annotation.gff3.gz \
    --refcds $ref_dir/dndscv/RefCDS_human_GRCh38_GencodeV18_recommended.rda \
    -N at31@sanger.ac.uk \
    -w ../work/ \
    -resume
)