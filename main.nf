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
pattern = file(params.pattern)
join_number = file(params.join_number)

include { ALIGN } from './subworkflow/align_bwa'

def getLibraryId( file ) {
    def filename = file.split(/\//)[-1]  // Extract filename from the full path
    def pattern = params.pattern  // Pattern to split the filename
    def join_number = params.join_number.toString().split(',').collect { it as int }  // Convert join_number to a list of integers

    // Split the filename by the given pattern
    def parts = filename.split(pattern)

    // Concatenate parts based on join_number indices
    def libraryId = join_number.collect { parts[it] }.join('_') 

    return libraryId
}


Channel
        .fromFilePairs( params.reads, flat: true )
        .map { prefix, file1, file2 -> tuple(getLibraryId(prefix), file1, file2) }
        .groupTuple()
        .set {reads}



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
	TARGETED_FAMILY(reads, ref_fa, ref_fai, known_snps_dbsnp_index, known_indels_index, known_snps_dbsnp, known_indels, target_bed)
}
