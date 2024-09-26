#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//include { ALIGN_READS } from "./process.nf"

def getLibraryId( file ) {
        file.split(/\//)[-1].split(/_/)[0]
}

Channel
        .fromFilePairs( params.reads, flat: true )
        .map { prefix, file1, file2 -> tuple(getLibraryId(prefix), file1, file2) }
        .groupTuple()
        .set {input}

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
	bwa-mem2.avx2 mem -t 20 -M -R '@RG\\tID:${samplename}\\tSM:${samplename}\\tPL:ILLUMINA' $params.ref $forward $reverse | gatk SortSam -I /dev/stdin -O ${samplename}.${params.timestamp}.sorted.bam --SORT_ORDER coordinate --CREATE_INDEX true

	cat <<-END_VERSIONS > versions.yml
	${task.process}\tbwa-mem2.avx2:\$(echo \$(bwa-mem2.avx2 version 2>&1) ); gatk:\$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
	END_VERSIONS
	"""
}

workflow {
    ALIGN_READS(input, ref_fa, ref_fai)
}
