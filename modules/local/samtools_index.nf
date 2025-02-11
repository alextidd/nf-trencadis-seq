// check bam is not truncated before proceeding + index
process samtools_index {
  tag "${meta.id}"
  label 'normal'

  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path(bam), path("${meta.id}.*am.*ai", arity: '1')

  script:
  """
  module load samtools-1.19/python-3.12.0 
  samtools quickcheck ${bam}
  samtools index -@ ${task.cpus} ${bam}
  """
}