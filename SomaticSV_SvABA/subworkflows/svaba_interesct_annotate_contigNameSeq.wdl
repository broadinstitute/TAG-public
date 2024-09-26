# svaba_intersect_annotate_contigNameSeq.wdl
# Author: Mohammad Tanhaemami
# Dr. Brian Crompton Lab
# Department of Pediatric Oncology
# Dana-Farber Cancer Institute

# This script annotates the vcf files that are the result of structural variation and indel analysis by assembly (SVABA).

workflow svaba_intersect_annotate_contigNameSeq {
	String SampleID
	File RefBed
	File InputVCF
	File InputContigBAM

	call intersectAnnotate {
	input:
		sampleID = SampleID,
		refBed = RefBed,
		inputVCF = InputVCF
	}
    call extractContigNameSeq {
	input:
		sampleID = SampleID,
		inputVCF = InputVCF,
		inputContigBam = InputContigBAM
	}
}

task intersectAnnotate {
	String sampleID
	File refBed
	File inputVCF

	Int? disk_size
    
    String vcf = basename(inputVCF, ".gz")
	String vcfBasename = basename(inputVCF)
	Boolean is_compressed = basename(inputVCF, ".vcf.gz") != vcfBasename

	command <<<
    
    # copy the input VCF file to the current working dir
    cp ${inputVCF} ${vcfBasename}
    
    # unzip the input vcf.gz file
	# gunzip -fc ${inputVCF} > ${sampleID}".svaba.unfiltered.somatic.sv.vcf"
    # decompress if the VCF file is compressed
    if [ ${is_compressed} = true ] 
    then
      gunzip -f ${vcfBasename}
    fi
	
	# find intersections between the sample file and the reference bed file
	bedtools intersect -header -wao -a ${vcf} -b ${refBed} > bufferIntersections.vcf
    
	# trim and remove the unwanted columns (preserver only the gene symbol and gene region)
	# ---------------------------------------------------------------------------------------
    # The below command is temporarily commented out because the current reference BED files do not have gene region information (intron or exon)
    # cut -f-11,15,16 bufferIntersections.vcf > bufferIntersections_annotated.vcf
    # Instead of the above command, I temporarily add the below command to only take the gene symbols from the current ref BED file.
    # Once the ref BED files are fixed and include intron and exon info, I will remove this one and use the command above.
    cut -f-11,15 bufferIntersections.vcf > bufferIntersections.annotated.vcf
    
    # remove duplicate lines that were created because of intersecting the VCF input with the reference bed file
    uniq bufferIntersections.annotated.vcf > bufferIntersections.annotated.unique.vcf
    
	# Add new column headers to the VCF file
	# sed -e '144s/$/\tgene/; 144s/$/\tregion/' bufferIntersections_annotated.vcf > ${sampleID}".svaba.unfiltered.somatic.sv.annotated.vcf"
    sed -e '144s/$/\tgene/' bufferIntersections.annotated.unique.vcf > ${sampleID}".svaba.unfiltered.somatic.sv.annotatedRegions.vcf"
    	
	>>>
	output {
	File svabaAnnotatedRegionsVCF = "${sampleID}.svaba.unfiltered.somatic.sv.annotatedRegions.vcf"
	}
 	runtime {
 		docker: "us.gcr.io/cromptonlab-001-mt/broadinstitute-gatk" 
 		memory: "2 GB"
   	cpu: 1
   	disks: "local-disk " + select_first([disk_size, 20]) + " HDD"
 	}

}

task extractContigNameSeq {
	String sampleID
	File inputVCF
	File inputContigBam

	Int? disk_size
	String? line = ""
	String? dollar = "$"
    
    String vcf = basename(inputVCF, ".gz")
	String vcfBasename = basename(inputVCF)
	Boolean is_compressed = basename(inputVCF, ".vcf.gz") != vcfBasename

	command <<<
    
    # copy the input VCF file to the current working dir
    cp ${inputVCF} ${vcfBasename}
    
    # unzip the input vcf.gz file
	# gunzip -fc ${inputVCF} > ${sampleID}".svaba.unfiltered.somatic.sv.vcf"
    # decompress if the VCF file is compressed
    if [ ${is_compressed} = true ] 
    then
      gunzip -f ${vcfBasename}
    fi
	
	# extract the contig sequences
	#sed -e '1,144d; s/.*SCTG=\(.*\)\;SECONDARY.*/\1/; s/.*SCTG=\(.*\)\;SPAN.*/\1/' ${sampleID}".svaba.unfiltered.somatic.sv.vcf" > "${sampleID}"_contigNames.txt
	grep "EVDNC" ${vcf} | awk -F "SCTG=" '{split($2,contig,";"); print contig[1]}' > "${sampleID}"_contigNames.txt
    
        
	#while IFS= read -r line
	#do
    echo -e "contig\tstrand\tseq" > "${sampleID}"_contigSequences.txt
	samtools view -F 0x10 ${inputContigBam} | cut -f 1,10 | sort -u | awk 'BEGIN{FS="\t"; OFS="\t"}{print $1,"+",$2}' >> "${sampleID}"_contigSequences.txt
    samtools view -f 0x10 ${inputContigBam} | cut -f 1,10 | sort -u | awk 'BEGIN{FS="\t"; OFS="\t"}{print $1,"-",$2}' >> "${sampleID}"_contigSequences.txt
	#done < "${sampleID}"_contigNames.txt
	
	>>>
	output {
	File contigNames = "${sampleID}_contigNames.txt"
	File contigSequences = "${sampleID}_contigSequences.txt"
	}
 	runtime {
 		docker: "us.gcr.io/cromptonlab-001-mt/broadinstitute-gatk" 
 		memory: "2 GB"
   	cpu: 1
   	disks: "local-disk " + select_first([disk_size, 20]) + " HDD"
 	}

}
