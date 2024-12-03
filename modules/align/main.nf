process ALIGN_READS {
  container 'jxprismdocker/bwa_avx2_gatk:latest'
  input:
    tuple val(samplename), file(forward), file(reverse)
    file(ref_fa)
    file(ref_fai)

  output:
    tuple val(samplename), file("${samplename}.${params.timestamp}.sorted.bam"), file("${samplename}.${params.timestamp}.sorted.bai")
    path "versions.yml", emit: versions

  script:
  """
  bwa-mem2.avx2 mem -t 20 -M -R '@RG\\tID:${samplename}\\tSM:${samplename}\\tPL:ILLUMINA' ${ref_fa} $forward $reverse | gatk SortSam -I /dev/stdin -O ${samplename}.${params.timestamp}.sorted.bam --SORT_ORDER coordinate --CREATE_INDEX true

  cat <<-END_VERSIONS > versions.yml
  ${task.process}\tbwa-mem2.avx2:\$(echo \$(bwa-mem2.avx2 version 2>&1) ); gatk:\$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
  END_VERSIONS
  """
}
