process BASE_RECALIBRATOR_WES{
  container 'jxprismdocker/prism_bwa_gatk'
  publishDir "$params.publishdir/$samplename", mode: 'copy', exclude: '*.yml'

  input:
    tuple val(samplename), file(sortedbam), file(sortedbam_index)
    file(ref_fa)
    file(ref_fai)
    file(known_snps_dbsnp_index)
    file(known_indels_index)
    file(known_snps_dbsnp)
    file(known_indels)
    file(target_bed)
  output:
    tuple val(samplename), file("${sortedbam.simpleName}.${params.timestamp}.recal_data.table")
    path "versions.yml", emit: versions

  script:
  """
  gatk BaseRecalibrator -R $params.ref -I $sortedbam --known-sites $known_snps_dbsnp --known-sites $known_indels -L $target_bed -O ${sortedbam.simpleName}.${params.timestamp}.recal_data.table
  
  cat <<-END_VERSIONS > versions.yml
  ${task.process}\tgatk:\$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
  END_VERSIONS
  """
}
