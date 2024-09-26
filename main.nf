#!/usr/bin/env nextflow
nextflow.enable.dsl=2

ref_fa = file(params.ref)
ref_fai = file(params.ref_fai)
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

workflow {
	MERGE_FASTQ(reads)
	//ALIGN_READS(MERGE_FASTQ.out, ref_fa, ref_fai)
}
