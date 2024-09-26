#!/usr/bin/env nextflow
nextflow.enable.dsl=2

ref_fa = file(params.ref)
ref_fai = file(params.ref_fai)
known_snps_dbsnp=file(params.known_snps_dbsnp)
known_snps_dbsnp_index=file(params.known_snps_dbsnp + '.tbi')
known_indels=file(params.known_indels)
known_indels_index=file(params.known_indels + '.tbi')
target_bed=file(params.target_bed)
target_bed_covered=file(params.target_bed_covered)
timestamp = file(params.timestamp)

include { MERGE_FASTQ } from './modules/merge_fastq'
include { ALIGN_READS } from './modules/align'

def getLibraryId( file ) {
        file.split(/\//)[-1].split(/_/)[0]
}

Channel
        .fromFilePairs( params.reads, flat: true )
        .map { prefix, file1, file2 -> tuple(getLibraryId(prefix), file1, file2) }
        .groupTuple()
        .set {reads}

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

/*workflow PRE-PROCESSING {
	MARK_DUPLICATES(ALIGN_READS.out[0])
        BASE_RECALIBRATOR_WES(MARK_DUPLICATES.out[0], ref_fa, ref_fai, known_snps_dbsnp_index, known_indels_index, known_snps_dbsnp, known_indels, target_bed)
        APPLY_BQSR_WES(MARK_DUPLICATES.out[0].join(BASE_RECALIBRATOR_WES.out[0]), ref_fa, ref_fai)
}

workflow VARIANT_CALLING {
	HAPLOTYPECALLER_WES(APPLY_BQSR_WES.out[0], ref_fa, ref_fai, known_snps_dbsnp_index, known_indels_index, known_snps_dbsnp, known_indels, target_bed)
        GENOTYPEGVCFS_WES(HAPLOTYPECALLER_WES.out[1].collect(),HAPLOTYPECALLER_WES.out[2].collect(), ref_fa, target_bed, params.proband)
	HARDFILTER_VARIANTS(GENOTYPEGVCFS_WES.out[0], ref_fa, ref_fai)
        DECOMPOSE_AND_NORMALIZE(HARDFILTER_VARIANTS.out[0], ref_fa, ref_fai)
}*/

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

workflow {
	TARGETED_FAMILY(reads, ref_fa, ref_fai, known_snps_dbsnp_index, known_indels_index, known_snps_dbsnp, known_indels, target_bed, ch_versions)
}
