process MERGE_FASTQ {
    input:
      tuple val(samplename), file(forward), file(reverse)
      
    output:
      tuple val(samplename), file("${samplename}_mergefastq_R1.fq.gz"), file("${samplename}_mergefastq_R2.fq.gz")

    script:
    """
    myvar=(${forward})
    myvar2=(${reverse})
    if ((\${#myvar[@]} > 1))
    then
      echo "Multiple lanes detected. Concatenating"
      cat ${forward} > ${samplename}_mergefastq_R1.fq.gz
      cat ${reverse} > ${samplename}_mergefastq_R2.fq.gz
    else
      echo "Single lane detected. Renaming symlink"
      mv \${myvar[0]} ${samplename}_mergefastq_R1.fq.gz
      mv \${myvar2[0]} ${samplename}_mergefastq_R2.fq.gz
    fi
    """
}
