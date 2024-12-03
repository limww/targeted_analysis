#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { MARK_DUPLICATES } from '../../modules/mark_duplicates'
include { BASE_RECALIBRATOR_WES } from '../../modules/base_recalibrator_wes'
include { APPLY_BQSR_WES } from '../../modules/apply_bqsr_wes'

workflow PREPROCESSING {
  take:
    ch_aligned_bam
    ref
    ref_fai
    known_snps_dbsnp_index
    known_indels_index
    known_snps_dbsnp
    known_indels
    target_bed
  main:
     ch_versions = Channel.empty() 
     MARK_DUPLICATES(ch_aligned_bam)
     ch_versions = ch_versions.mix(MARK_DUPLICATES.out.versions)
     BASE_RECALIBRATOR_WES(MARK_DUPLICATES.out[0], ref_fa, ref_fai, known_snps_dbsnp_index, known_indels_index, known_snps_dbsnp, known_indels, target_bed)
     ch_versions = ch_versions.mix(BASE_RECALIBRATOR_WES.out.versions)
     APPLY_BQSR_WES(MARK_DUPLICATES.out[0].join(BASE_RECALIBRATOR_WES.out[0]), ref_fa, ref_fai)
     ch_versions = ch_versions.mix(APPLY_BQSR_WES.out.versions)
  emit:
    mark_dup_bam = MARK_DUPLICATES.out[0]
    bqsr_recal_table = BASE_RECALIBRATOR_WES.out[0]
    bqsr_bam = APPLY_BQSR_WES.out[0]
    versions = ch_versions
}
