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
               germline_sv_bed = germline_sv_bed,
               xomere_bed = xomere_bed,
               target_bed = target_bed
    }
    output {
        Svaba.*
    }
}