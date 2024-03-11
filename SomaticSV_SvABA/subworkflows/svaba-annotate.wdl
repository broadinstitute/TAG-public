# svaba-annotate-Walaj.wdl
# Author: Mohammad Tanhaemami
# Dr. Brian Crompton Lab
# Department of Pediatric Oncology
# Dana-Farber Cancer Institute

# This script annotates input vcf files that are the outputs of structural variation and indel analysis by assembly (SVABA) workflow.

workflow SvabaAnnotationWalajRscript {
	String SampleID
	File InputVCF
    String? Style
    String GenomeBuild
    String? GeneDataBase
    File Svaba_Annotate_Rscript

	call AnnotateSvabaWalajR {
	input:
		sampleID = SampleID,
		inputVCF = InputVCF,
        style = Style,
        genomeBuild = GenomeBuild ,
        geneDataBase = GeneDataBase,
        svaba_annotate_Rscript = Svaba_Annotate_Rscript
	}
}

task AnnotateSvabaWalajR {
	String sampleID
	File inputVCF
    String? style
    String genomeBuild
    String? geneDataBase
    File svaba_annotate_Rscript

	Int? disk_size

	command <<<
    	Rscript ${svaba_annotate_Rscript} --input ${inputVCF} --genome ${genomeBuild} --output ${sampleID}
	>>>
	output {
	 File annotatedVCFwalajRscript = "${sampleID}.svaba.unfiltered.somatic.sv.annotated.vcf"
     File full_table_r_object = "${sampleID}.full_vcf_object.rds"
	}
 	runtime {
 		#docker: "us-central1-docker.pkg.dev/cromptonlab-001-mt/cromptonlab-docker-repo/jnktsj-svaba-v0.2.1:1"
        #docker: "us.gcr.io/broad-dsp-gcr-public/terra-jupyter-r"
        docker: "us-central1-docker.pkg.dev/cromptonlab-001-mt/cromptonlab-docker-repo/bioconductor-with-additional-r-packages:1"
 		memory: "8 GB"
   	cpu: 1
   	disks: "local-disk " + select_first([disk_size, 20]) + " HDD"
 	}

}
