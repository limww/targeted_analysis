include { MERGE_FASTQ } from '../../modules/merge_fastq'
include { ALIGN_READS } from '../../modules/align'

workflow ALIGN {
	take:
		reads
		ref
		ref_fai
	main:
		ch_versions = Channel.empty()
		MERGE_FASTQ(reads)
		ALIGN_READS(MERGE_FASTQ.out, ref_fa, ref_fai)
	emit: 
		aligned_bam = ALIGN_READS.out[0]
		ch_versions = ch_versions.mix(ALIGN_READS.out.versions)
}
