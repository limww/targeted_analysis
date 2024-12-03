#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { MERGE_FASTQ } from '../../modules/merge_fastq'
include { ALIGN_READS } from '../../modules/align'

workflow BWA_ALIGN {
	take:
		reads
		ref_fa
		ref_fai
	main:
		ch_versions = Channel.empty()
		MERGE_FASTQ(reads)
		ALIGN_READS(MERGE_FASTQ.out, ref_fa, ref_fai)
		ch_versions = ch_versions.mix(ALIGN_READS.out.versions)
	emit: 
		aligned_bam = ALIGN_READS.out[0]
		versions = ch_versions
		
}
