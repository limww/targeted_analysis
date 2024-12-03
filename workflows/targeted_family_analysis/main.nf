#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { ALIGN } from "../../subworkflow/align_bwa"
  
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
		ALIGN(reads, ref, ref_fai)
		//PRE-PROCESSING(ref_fa, ref_fai, known_snps_dbsnp_index, known_indels_index, known_snps_dbsnp, known_indels, target_bed)
		//VARIANT_CALLING(ref_fa, ref_fai, known_snps_dbsnp_index, known_indels_index, known_snps_dbsnp, known_indels, target_bed, params.proband)
	emit:
		ALIGN.out.aligned_bam
		ALIGN.out.ch_versions
}
