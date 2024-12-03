#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { BWA_ALIGN } from '../../subworkflow/align_bwa'
INCLUDE { GATK_PRACTICES }  from '../../subworkflow/gatk_practices'
workflow TARGETED_FAMILY{
	take: 
		reads
		ref
		ref_fai
		known_snps_dbsnp_index
		known_indels_index
		known_snps_dbsnp
		known_indels
		target_bed
		ch_versions
	main:
		ch_versions = Channel.empty()
		BWA_ALIGN(reads, ref, ref_fai)
		ch_versions = ch_versions.mix(BWA_ALIGN.out.versions)
		
		ch_aligned_bam = BWA_ALIGN.out.aligned_bam
		GATK_PRACTICES(ch_aligned_bam, ref_fa, ref_fai, known_snps_dbsnp_index, known_indels_index, known_snps_dbsnp, known_indels, target_bed)
		ch_versions = ch_versions.mix(GATK_PRACTICES.out.versions)

		//VARIANT_CALLING(ref_fa, ref_fai, known_snps_dbsnp_index, known_indels_index, known_snps_dbsnp, known_indels, target_bed, params.proband)
	emit:
		BWA_ALIGN.out.aligned_bam
		GATK_PRACTICES.bqsr_recal_table
		GATK_PRACTICES.bqsr_bam
		versions = ch_versions
}
