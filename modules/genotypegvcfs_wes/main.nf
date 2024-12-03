process GENOTYPEGVCFS_WES{
  container 'jxprismdocker/prism_bwa_gatk'
  publishDir "$params.publishdir/jointcalling", mode: 'copy', exclude: '*.yml'
  
  input:
    file(gvcf)
    file(gvcf_tbi)
    file(ref_fai)
    file(target_bed)
    stdin samplename
  
  output:
    tuple val(samplename), file("${samplename}.${params.timestamp}.raw.vcf.gz"), file("${samplename}.${params.timestamp}.raw.vcf.gz.tbi")
    path "versions.yml", emit: versions
  
  script:
  samplename = params.proband
  """
  echo *.g.vcf.gz | tr " " "\n" > gvcf.list
  gatk CombineGVCFs -R $params.ref_fa -V gvcf.list -O ${samplename}.combinedgvcf.vcf.gz
  gatk GenotypeGVCFs -R $params.ref_fa -V ${samplename}.combinedgvcf.vcf.gz  -O ${samplename}.${params.timestamp}.raw.vcf.gz -L $target_bed
  
  cat <<-END_VERSIONS > versions.yml
  ${task.process}\tgatk:\$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
  END_VERSIONS
  """
}
