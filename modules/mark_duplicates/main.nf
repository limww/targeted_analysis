process MARK_DUPLICATES{
  container 'jxprismdocker/prism_bwa_gatk'
  input:
    tuple val(samplename), file(bamfile), file(bamindex)

  output:
    tuple val(samplename), file("${bamfile.simpleName}.${params.timestamp}.sorted.rmdup.bam"), file("${bamindex.simpleName}.${params.timestamp}.sorted.rmdup.bai")
    tuple val(samplename), file("${bamfile.simpleName}.${params.timestamp}.marked_dup_metrics.txt")
	  path "versions.yml", emit: versions

  script:
  """
  gatk MarkDuplicates -I $bamfile -O ${bamfile.simpleName}.${params.timestamp}.sorted.rmdup.bam --CREATE_INDEX true --METRICS_FILE ${bamfile.simpleName}.${params.timestamp}.marked_dup_metrics.txt

  cat <<-END_VERSIONS > versions.yml
  ${task.process}\tgatk:\$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
  END_VERSIONS
  """
}
