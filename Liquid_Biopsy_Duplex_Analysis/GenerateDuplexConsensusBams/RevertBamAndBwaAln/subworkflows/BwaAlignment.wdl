workflow BwaAlignmentTest {
	call BwaAlignment
}

task BwaAlignment {
	File refFasta
	File refFastaIndex
	File refFastaDict
    File ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
	File firstEndFastq
    String fq1 = basename(firstEndFastq)
    String basename1 = basename(firstEndFastq, ".fastq.gz")
	File secondEndFastq
    String fq2 = basename(secondEndFastq)
    String basename2 = basename(secondEndFastq, ".fastq.gz")
	String sampleName
    String gitc_docker
	Int memoryGb
	Int diskSpaceGb
	Int cpu

	command <<<
    
    	mv ${firstEndFastq} ./${fq1}
        mv ${secondEndFastq} ./${fq2}
        
        /usr/gitc/bwa aln -q 5 -l 32 -k 2 -t ${cpu} -o 1 ${refFasta} ./${fq1} -f ./${basename1}.sai
        export bwa_cmd="/usr/gitc/bwa aln -q 5 -l 32 -k 2 -t "${cpu}" -o 1 "${refFasta}" ./"${fq1}" -f ./"${basename1}".sai\;"
        
		/usr/gitc/bwa aln -q 5 -l 32 -k 2 -t ${cpu} -o 1 ${refFasta} ./${fq2} -f ./${basename2}.sai
        export bwa_cmd=$bwa_cmd" /usr/gitc/bwa aln -q 5 -l 32 -k 2 -t "${cpu}" -o 1 "${refFasta}" ./"${fq2}" -f ./"${basename2}".sai\;"

		/usr/gitc/bwa sampe -P ${refFasta} ./${basename1}.sai ./${basename2}.sai ./${fq1} ./${fq2} -f ./${sampleName}.aligned.sam
        export bwa_cmd=$bwa_cmd" /usr/gitc/bwa sampe -P "${refFasta}" ./"${basename1}".sai ./"${basename2}".sai ./"${fq1}" ./"${fq2}" -f ./"${sampleName}".aligned.sam"
        echo $bwa_cmd > bwa_cmd.txt
        
        samtools sort -n ${sampleName}.aligned.sam -o ${sampleName}.aligned.bam
        
	>>>

	output {
        File raw_aligned_bam = "${sampleName}.aligned.bam"
        String bwa_command = read_string("bwa_cmd.txt")
	}

	runtime {
		docker: gitc_docker
		memory: "${memoryGb} GB"
		cpu: "${cpu}"
		disks: "local-disk ${diskSpaceGb} HDD"
	}

}