import "./subworkflows/svaba_interesct_annotate_contigNameSeq.wdl" as SvABA_annotate_extractContigs
import "./subworkflows/snpEffCaller_extractFields.wdl" as snpEffCaller
import "./subworkflows/svaba-annotate.wdl" as walaj_annotate_SvABA
import "./subworkflows/MergeSvabaAnnotations.wdl" as mergeSvabaAnnotations

task Svaba {
    Int disk_size
    Int cpu_num
    Int memory_size

    # sample inputs
    String sample_name
    File tumor_bam
    File tumor_bam_index
    File control_bam
    File control_bam_index
    
    # reference files
    File reference_fasta
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa        

    File dbsnp_vcf
    File germline_sv_bed
    File xomere_bed
    String? svaba_extra_args

    File? target_bed

    command <<<
        # sample links
        ln -vs ${tumor_bam} ${sample_name}_tumor.bam
        ln -vs ${tumor_bam_index} ${sample_name}_tumor.bai
        ln -vs ${control_bam} ${sample_name}_control.bam
        ln -vs ${control_bam_index} ${sample_name}_control.bai

        # reference links
        ln -vs ${reference_fasta} reference.fasta
        ln -vs ${reference_fai} reference.fasta.fai
        ln -vs ${reference_amb} reference.fasta.amb
        ln -vs ${reference_ann} reference.fasta.ann
        ln -vs ${reference_bwt} reference.fasta.bwt
        ln -vs ${reference_pac} reference.fasta.pac
        ln -vs ${reference_sa} reference.fasta.sa

        /svaba/src/svaba/svaba run -t ${tumor_bam} \
                                   -n ${control_bam} \
                                   -G reference.fasta \
                                   -p ${cpu_num} \
                                   -D ${dbsnp_vcf} \
                                   -V ${germline_sv_bed} \
                                   -B ${xomere_bed} \
                                   ${"-k " + target_bed} \
                                   ${svaba_extra_args} \
                                   -a ${sample_name} &&
        
        # bgzip unfiltered VCFs
        bgzip ${sample_name}.svaba.unfiltered.somatic.sv.vcf
        bgzip ${sample_name}.svaba.unfiltered.somatic.indel.vcf
        bgzip ${sample_name}.svaba.unfiltered.germline.sv.vcf
        bgzip ${sample_name}.svaba.unfiltered.germline.indel.vcf
        
        # extract the number of detected SV events
        if [[ `grep -v "#" ${sample_name}.svaba.somatic.sv.vcf | wc -l` != 0 ]]
        then
            grep -v "#" ${sample_name}.svaba.somatic.sv.vcf | \
            awk '{ split($5, tt, /\[|\]/); split(tt[2], chr, ":");
                   split($3, ID, ":");
                   if($1 == chr[1]) {
                      print ID[1]"\tIntra"
                   } else {
                      print ID[1]"\tInter"
                   } }' | sort -u > called_sv_type
            cat called_sv_type | wc -l > somatic_sv_count
            awk '{ if($2 == "Intra"){print} }' called_sv_type | wc -l > intrachromosomal
            awk '{ if($2 == "Inter"){print} }' called_sv_type | wc -l > interchromosomal
        else
            echo "0" > somatic_sv_count
            echo "0" > intrachromosomal
            echo "0" > interchromosomal
        fi
    >>>
    runtime {
        docker: "jnktsj/svaba-v0.2.1:1"
        memory: "${memory_size} GB"
        cpu: "${cpu_num}"
        disks: "local-disk ${disk_size} HDD"
        preemptible: 0
    }
    output {
        # VCF outputs
        File somatic_sv = "${sample_name}.svaba.somatic.sv.vcf"
        File somatic_indel = "${sample_name}.svaba.somatic.indel.vcf"
        File germline_sv = "${sample_name}.svaba.germline.sv.vcf"
        File germline_indel = "${sample_name}.svaba.germline.indel.vcf"

        # Unfiltered VCF outputs
        File somatic_unfiltered_sv = "${sample_name}.svaba.unfiltered.somatic.sv.vcf.gz"
        File somatic_unfiltered_indel = "${sample_name}.svaba.unfiltered.somatic.indel.vcf.gz"        
        File germline_unfiltered_sv = "${sample_name}.svaba.unfiltered.germline.sv.vcf.gz"
        File germline_unfiltered_indel = "${sample_name}.svaba.unfiltered.germline.indel.vcf.gz"

        # Other outputs
        File discordant_read = "${sample_name}.discordant.txt.gz"
        File contig_bam = "${sample_name}.contigs.bam"
        File contig_alignment = "${sample_name}.alignments.txt.gz"
        File raw_unfiltered_variant = "${sample_name}.bps.txt.gz"
        File svaba_log = "${sample_name}.log"

        # SV stats
        String somatic_sv_count = read_string("somatic_sv_count")
        String intrachromosomal = read_string("intrachromosomal")
        String interchromosomal = read_string("interchromosomal")
    }
}


task IntervalToBed {
    File? target_interval
    String output_bed = "target_interval.bed"

    command <<<
        if [[ "${target_interval}" == *.interval_list ]]
        then
            grep -v "^@" ${target_interval} | \
            awk '{if($3>$2){$2-=1; print $1,$2,$3,$5}}' OFS="\t" > ${output_bed}
        else
            cp ${target_interval} ${output_bed}
        fi
    >>>

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.2.5-1486412288"
        disks: "local-disk 1 HDD"
        memory: "1 GB"
    }
    output {
        File target_bed = "${output_bed}"
    }
}
								  

workflow SomaticRearrangement {

    Int disk_size
    Int memory_size
    Int cpu_num

    String sample_name
    File tumor_bam
    File tumor_bam_index
    File control_bam
    File control_bam_index
    
    File reference_fasta
    File reference_fai
    File reference_amb
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_sa

    File dbsnp_vcf
    File germline_sv_bed
    File xomere_bed
    File? target_interval
    String? svaba_extra_args
    
    # Svaba annotation inputs
    File annotation_reference_bed
    File snpEff_jar
    File snpSift_jar
    File configFile
    File databaseFile
    String? Style
    String GenomeBuild
    String? GeneDataBase
    File Svaba_Annotate_Rscript
	File Merge_Annotations_Script
    
    Boolean annotate_contigs = true
    Boolean extract_contig_sequences = true
    Boolean merge_annotations = true

    if (defined(target_interval)) {
        call IntervalToBed {
            input: target_interval = target_interval
        }
    }

    File? target_bed = IntervalToBed.target_bed

    call Svaba {
        input: disk_size = disk_size,
               memory_size = memory_size,
               cpu_num = cpu_num,
               sample_name = sample_name,
               tumor_bam = tumor_bam,
               tumor_bam_index = tumor_bam_index,
               control_bam = control_bam,
               control_bam_index = control_bam_index,
               reference_fasta = reference_fasta,
               reference_fai = reference_fai,
               reference_amb = reference_amb,
               reference_ann = reference_ann,
               reference_bwt = reference_bwt,
               reference_pac = reference_pac,
               reference_sa = reference_sa,
               dbsnp_vcf = dbsnp_vcf,
               svaba_extra_args = svaba_extra_args,
               germline_sv_bed = germline_sv_bed,
               xomere_bed = xomere_bed,
               target_bed = target_bed
    }
    if(annotate_contigs){
    	call SvABA_annotate_extractContigs.intersectAnnotate as SvabaAnnotation {
    		input:
            	sampleID = sample_name,
        	   	refBed = annotation_reference_bed,
                inputVCF = Svaba.somatic_unfiltered_sv
    	}
        call snpEffCaller.snpEffCaller as SvabaAnnotationSnpEff {
        	input:
            	snpEff_jar = snpEff_jar,
                snpSift_jar = snpSift_jar,
           		inputVCF = Svaba.somatic_unfiltered_sv,
           		configFile = configFile,
           		databaseFile = databaseFile,
           		outputName = sample_name
        }
        call walaj_annotate_SvABA.AnnotateSvabaWalajR as SvabAnnotationWalaj {
        	input:
            	sampleID = sample_name,
				inputVCF = Svaba.somatic_unfiltered_sv,
        		style = Style,
        		genomeBuild = GenomeBuild ,
        		geneDataBase = GeneDataBase,
        		svaba_annotate_Rscript = Svaba_Annotate_Rscript
        }
    }
    if(extract_contig_sequences){
    	call SvABA_annotate_extractContigs.extractContigNameSeq as SvabaContigExtraction {
        input:
        	sampleID = sample_name,
        	inputVCF = Svaba.somatic_unfiltered_sv,
            inputContigBam = Svaba.contig_bam
        }
    }
    if (merge_annotations){
    	call mergeSvabaAnnotations.mergeAnnotationsTask as mergeAnnotations {
        input:
        	output_basename = sample_name,
            merge_annotations_script = Merge_Annotations_Script,
            svaba_vcf = Svaba.somatic_unfiltered_sv,
            snpsift_table = SvabaAnnotationSnpEff.vcfOutExtractFields,
            annot_table = SvabAnnotationWalaj.annotatedVCFwalajRscript,
            contig_seqs = SvabaContigExtraction.contigSequences
        }
    }
    
    output {
        Svaba.*
        File? annotated_sv_vcf = SvabaAnnotation.svabaAnnotatedRegionsVCF
        File? annotated_sv_vcf_snpEff = SvabaAnnotationSnpEff.vcfOut
        File? annotated_sv_vcf_snpEff_extractedFields = SvabaAnnotationSnpEff.vcfOutExtractFields
        File? annotated_sv_vcf_walaj = SvabAnnotationWalaj.annotatedVCFwalajRscript
        File? walaj_r_object_table = SvabAnnotationWalaj.full_table_r_object
        File? contig_names = SvabaContigExtraction.contigNames
        File? contig_sequences = SvabaContigExtraction.contigSequences
        File? final_annotation = mergeAnnotations.final_annotation
    }
}