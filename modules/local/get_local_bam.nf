// symlink local bam
process get_local_bam {
  tag "${meta.id}"
  maxForks 10
  label 'normal4core'

  input:
  tuple val(meta), val(bam)

  output:
  tuple val(meta), path("${meta.id}.*am", arity: '1')

  script:
  def bam_ext = bam.getName().substring(bam.getName().lastIndexOf('.') + 1)
  """
  ln -s ${bam} ${meta.id}.${bam_ext}
  """
}