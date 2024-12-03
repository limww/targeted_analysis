process APPLY_BQSR_WES {
  container 'jxprismdocker/prism_bwa_gatk'
  publishDir "$params.publishdir/$samplename", mode: 'copy', exclude: '*.yml'
  input:
    tuple val(samplename), file(sortedbam), file(sortedbam_index), file(recal_table)
    file(ref_fa)
    file(ref_fai)
  output:
    tuple val(samplename), file("${sortedbam.simpleName}.${params.timestamp}.BQSR.bam"), file("${sortedbam.simpleName}.${params.timestamp}.BQSR.bai")
    path "versions.yml", emit: versions
  
  script:
  """
  gatk ApplyBQSR -R $params.ref_fa -I $sortedbam --bqsr-recal-file ${recal_table} -O ${sortedbam.simpleName}.${params.timestamp}.BQSR.bam

  cat <<-END_VERSIONS > versions.yml
  ${task.process}\tgatk:\$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
  END_VERSIONS
  """
}
