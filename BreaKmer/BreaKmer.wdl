# This workflow is originally from CCGD DFCI and modified as a Terra/FireCloud version by TAG
import 
"https://raw.githubusercontent.com/broadinstitute/TAG-public/99c8c2ec6d7d6b798b62435ed260b0a6decdacd2/BreaKmer/breakmer_subworkflows/AlignAndHardClipBam.wdl" 
as realign_hardclip
import "https://github.com/broadinstitute/TAG-public/blob/99c8c2ec6d7d6b798b62435ed260b0a6decdacd2/BreaKmer/breakmer_subworkflows/AnnotateBreaKmer.wdl" as 
AnnotateBreaKmer


workflow BreakmerAnalysis{
    File ref_fasta
    File ref_fasta_index
    File ref_2bit

    # These files are target specific reference files.
    # Once you have run BreaKmer with a target list, the workflow
    # should generate a tar gz file for the targets.
    File? ref_kmer_tar

    # BED file that contains repeat regions in the genome:
    #    [1] chromosome
    #    [2] start
    #    [3] end
    #    [4] repeat name
    File repeat_mask_bed

    # Gene annotation in the refGene format, see "Gene Preditions (Extended)"
    # for more details: https://genome.ucsc.edu/FAQ/FAQformat.html#format9
    # BreaKmer reads:
    #   3rd column for choromosome names
    #   5th column for start
    #   6th column for end
    #  13th column for gene name
    File refgene_annotation

    # Other required and optional parameters for BreaKmer
    File targets_bed
    File cutadapt_config
    File? gene_list
    File? other_regions_bed

    # Default parameters for k-mer and indel sizes are set to 15
    Int? kmer_size
    Int? indel_size

    # Override paths
    String? docker
    String? cutadapt_path
    String? blat_path
    String? jellyfish_path

    # Scatter number for BreaKmer runs
    Int? scatter_num

    # Run MarkDuplicates
    Boolean? run_markduplicates
    Boolean? remove_duplicates_from_bam
    Boolean remove_duplicates = select_first([remove_duplicates_from_bam, true])
    
    # Run Downsampling
    Float? downsample_frac
    File ds_stats_python
    Boolean? run_downsample_bam
    Boolean run_ds_bam = select_first([run_downsample_bam, true])
    
    # Run Realignment with HardClipping
    Boolean? run_realign_hardclipping
    Boolean run_realign_hc = select_first([run_realign_hardclipping, true])
    Boolean fix_mate = if (run_ds_bam || run_md) then true else false
    
    Boolean run_md = if (run_ds_bam || run_realign_hc) then false else select_first([run_markduplicates, true])
    
    # Annotate BreaKmer Results
    Boolean? annotate_breakmer_sv = true
    Boolean run_annotate_breakmer = select_first([annotate_breakmer_sv, true])
    String? annotate_breakmer_docker
    File annotation_bed

    String? picard_docker
    File? picard_override

    String sample_name
    File input_bam
    File input_bam_index

    # Additional disk size to pad, default value is 10 GB
    Int? disk_pad
    Int ref_2bit_size = if defined(ref_2bit) then ceil(size(ref_2bit, "GB")) else 0
    Int ref_kmer_size = if defined(ref_kmer_tar) then ceil(size(ref_kmer_tar, "GB")) else 0
    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + ref_2bit_size)
    Int original_bam_size = ceil(size(input_bam, "GB") + size(input_bam_index, "GB"))
    Int aux_size = ceil(size(targets_bed, "GB")/select_first([scatter_num, 20]) + size(repeat_mask_bed, "GB") + size(cutadapt_config, "GB"))

    if(run_md) {
        call MarkDuplicates {
            input:
                input_bam = input_bam,
                input_bam_index = input_bam_index,
                remove_duplicates_from_bam = remove_duplicates,
                picard_docker = picard_docker,
                picard_override = picard_override,
                disk_size = original_bam_size * 3 + select_first([disk_pad, 10])
        }
    }

    if(run_ds_bam) {
        Int input_bam_size = ceil(size(input_bam, "GB") + size(input_bam_index, "GB"))

        call DynamicDownsampleBam {
            input:
                input_bam = input_bam,
                input_bam_index = input_bam_index,
                downsample_frac = downsample_frac,
                ds_stats_python = ds_stats_python,
                targets_bed = targets_bed,
                picard_docker = picard_docker,
                picard_override = picard_override,
                disk_size = input_bam_size * 3 + select_first([disk_pad, 10])
        }
    }
    
    if(run_realign_hc){
    	call realign_hardclip.AlignAndHardClipBam as realign_hardclip{
        	input: input_bam = select_first([DynamicDownsampleBam.output_bam, input_bam]),
            	   base_file_name = sample_name,
                   remove_duplicates = remove_duplicates,
                   fix_mate = fix_mate
        }
    }

    File bam = select_first([realign_hardclip.output_bam, DynamicDownsampleBam.output_bam, MarkDuplicates.output_bam, input_bam])
    File bam_index = select_first([realign_hardclip.output_bam_index, DynamicDownsampleBam.output_bam_index, MarkDuplicates.output_bam_index, input_bam_index])
    Int bam_size = ceil(size(bam, "GB") + size(bam_index, "GB"))

    # Default setting splits a target file into 20 small files
    call SplitTargets {
        input:
            targets_bed = targets_bed,
            scatter_num = select_first([scatter_num, 20]),
            disk_size = ceil(size(targets_bed, "GB")) + 2
    }


    scatter (subtarget_file in SplitTargets.target_files) {
        call RunBreakmer{
            input:
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_2bit = ref_2bit,
                ref_kmer_tar = ref_kmer_tar,
                repeat_mask_bed = repeat_mask_bed,
                refgene_annotation = refgene_annotation,
                targets_bed = subtarget_file,
                cutadapt_config = cutadapt_config,
                gene_list = gene_list,
                other_regions_bed = other_regions_bed,
                kmer_size = select_first([kmer_size, 15]),
                indel_size = select_first([indel_size, 15]),
                docker = select_first([docker, "us.gcr.io/tag-team-160914/breakmer:0.0.6"]),
                cutadapt_path = select_first([cutadapt_path, "/usr/local/bin/cutadapt"]),
                blat_path = select_first([blat_path, "/apps/blat_35x1/blat"]),
                jellyfish_path = select_first([jellyfish_path, "/apps/jellyfish-1.1.11/bin/jellyfish"]),
                sample_name = sample_name,
                input_bam = bam,
                input_bam_index = bam_index,
                disk_size = ref_size + bam_size + aux_size + select_first([disk_pad, 10])
        }
    }

    call GatherBreakmerResults {
        input:
            sample_name = sample_name,
            summary_files = RunBreakmer.output_summary,
            sv_files = RunBreakmer.output_sv
    }
    if(run_annotate_breakmer){
    	call AnnotateBreaKmer.annotate_breakmer_sv as AnnotateBreaKmer {
    		input: 
            	annotation_bed = annotation_bed,
                breakmer_sv = GatherBreakmerResults.output_sv,
                docker_override = annotate_breakmer_docker,
                output_basename = sample_name                
    	}
	}
    
    output {
        File output_summary = GatherBreakmerResults.output_summary
        File output_sv = select_first([AnnotateBreaKmer.annotated_sv, GatherBreakmerResults.output_sv])
        Int total_sv_detected = GatherBreakmerResults.total_sv_detected

        File? deduplicated_bam = MarkDuplicates.output_bam
        File? deduplicated_bam_index = MarkDuplicates.output_bam_index
        File? duplicate_metrics =  MarkDuplicates.duplicate_metrics

        File? downsample_bam = DynamicDownsampleBam.output_bam
	File? downsample_bam_index = DynamicDownsampleBam.output_bam_index
        File? downsample_stats = DynamicDownsampleBam.downsample_stats
        Float? downsample_fraction = DynamicDownsampleBam.downsample_fraction
        
        File? realigned_bam = realign_hardclip.output_bam
        File? realigned_bam_index = realign_hardclip.output_bam_index
    }
}


task MarkDuplicates {
    File input_bam
    File input_bam_index
    Boolean remove_duplicates_from_bam
    String? picard_docker
    File? picard_override

    Int? preemptible_attempts
    Int? max_retries
    Int disk_size
    Int? memory

    Int mem = select_first([memory, 8])
    Int compute_mem = mem * 1000 - 500

    String out_basename = basename(input_bam, ".bam")

    command {
        set -e
        
        export PICARD_JAR=${default="/usr/gitc/picard.jar" picard_override}

        java -Xmx${compute_mem}m -jar $PICARD_JAR MarkDuplicates \
            INPUT=${input_bam} \
            OUTPUT=${out_basename}.markdup.bam \
            M=${out_basename}.duplicate_metrics \
            CREATE_INDEX=true ${true="REMOVE_DUPLICATES=true" false=" " remove_duplicates_from_bam}
    }
    output {
        File output_bam = "${out_basename}.markdup.bam"
        File output_bam_index = "${out_basename}.markdup.bai"
        File duplicate_metrics = "${out_basename}.duplicate_metrics"
    }
    runtime {
        docker: select_first([picard_docker, "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"])
        memory: mem + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        maxRetries: select_first([max_retries, 1])
    }
}


task DynamicDownsampleBam {
    File input_bam
    File input_bam_index
    String? picard_docker
    File? picard_override

    File targets_bed
    File ds_stats_python
    Float? downsample_frac
    Boolean run_dynamic_downsample = if defined(downsample_frac) then false else true
    Float ds_frac_default = select_first([downsample_frac, 1.0])

    Int? preemptible_attempts
    Int? max_retries
    Int disk_size
    Int? memory

    Int mem = select_first([memory, 8])
    Int compute_mem = mem * 1000 - 500

    String out_basename = basename(input_bam, ".bam")

    command <<<
        set -e

        DS_FRAC=${ds_frac_default}
        export PICARD_JAR=${default="/usr/gitc/picard.jar" picard_override}

        # write the default setting of downsampling parameters
        echo $DS_FRAC | awk '{print "[DOWNSAMPLE]\ndownsample_fraction:\t"$0}' > "${out_basename}.ds_stats.txt"

        if [ ${run_dynamic_downsample} = true ]
        then
            samtools idxstats ${input_bam} > idx_stats.out
            python ${ds_stats_python} ${targets_bed} idx_stats.out > "${out_basename}.ds_stats.txt"

            DS_FRAC=$(grep ^downsample_fraction "${out_basename}.ds_stats.txt" | cut -f2)
        fi

        if [ $DS_FRAC == 1.0 ]
        then
            cp ${input_bam} "${out_basename}.ds.bam"
            cp ${input_bam_index} "${out_basename}.ds.bai"
        else
            java -Xmx${compute_mem}m -jar $PICARD_JAR DownsampleSam \
                INPUT=${input_bam} \
                OUTPUT="${out_basename}.ds.bam" \
                P=$DS_FRAC \
                CREATE_INDEX=true
        fi

        echo $DS_FRAC > "ds_param.txt"
    >>>
    output {
        File output_bam = "${out_basename}.ds.bam"
	File output_bam_index = "${out_basename}.ds.bai"
        File downsample_stats = "${out_basename}.ds_stats.txt"
        Float downsample_fraction = read_float("ds_param.txt")
    }
    runtime {
        docker: select_first([picard_docker, "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"])
        memory: mem + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        maxRetries: select_first([max_retries, 1])
    }
}

task SplitTargets {
    File targets_bed
    Int scatter_num
    String out_basename = basename(targets_bed, ".bed")
    
    Int? preemptible_attempts
    Int? max_retries
    Int disk_size
    
    command <<<
        set -e -o pipefail
        
        # Count the number of targets
        TARGETS=`grep -c ^[A-Za-z0-9] ${targets_bed}`

        # Compute a number of targets per split file
        SPLIT_NUM=`echo $(( $(( $TARGETS / ${scatter_num} )) + $(( $TARGETS % ${scatter_num} > 0)) ))`
        
        sort -V ${targets_bed} > sorted.bed
        split -l $SPLIT_NUM sorted.bed "`basename ${targets_bed} .bed`."
    >>>
    output {
        Array[File] target_files = glob("${out_basename}.*")
    }
    runtime {
        docker: "ubuntu:16.04"
        memory: "4 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        maxRetries: select_first([max_retries, 1])    
    }
}


task RunBreakmer {
    File input_bam
    File input_bam_index

    File ref_fasta
    File ref_fasta_index
    File ref_2bit
    File? ref_kmer_tar
    File repeat_mask_bed
    File refgene_annotation
    File targets_bed
    File cutadapt_config
    File? gene_list
    File? other_regions_bed

    String cutadapt_path
    String blat_path
    String jellyfish_path
    String sample_name

    Int kmer_size
    Int indel_size

    # Runtime variables
    String docker
    Int? cpu
    Int? memory_gb
    Int? preemptible_attempts
    Int? max_retries
    Int disk_size

    String targets_basename = basename(targets_bed, ".bed")
    String config_file = sample_name + "_breakmer.config"
    String output_dir = sample_name + "_breakmer"
    String output_log = sample_name + "_breakmer.log"

    command {
        set -e -o pipefail

        # Prepare reference k-mer directory
        mkdir -p ${targets_basename}
        if [[ -f "${ref_kmer_tar}" ]]; then
            tar -xzvf ${ref_kmer_tar} --directory ${targets_basename} --strip-components=1
        fi
        REFDATA=`readlink -f ${targets_basename}`
        mkdir scratch
        SCRATCH=`readlink -f scratch`
        mkdir -p ${output_dir}
        OUTPUT=`readlink -f ${output_dir}`

        # Create configuration file
        echo "analysis_name=${sample_name}" > ${config_file}
        echo "sample_bam_file=${input_bam}" >> ${config_file}
        echo "analysis_dir=$OUTPUT" >> ${config_file}
        echo "targets_dir=$SCRATCH" >> ${config_file}
        echo "targets_bed_file=${targets_bed}" >> ${config_file}
        echo "gene_annotation_file=${refgene_annotation}" >> ${config_file}

        if [[ -f "${other_regions_bed}" ]]; then
            echo "other_regions_file=${other_regions_bed}" >> ${config_file}
        fi

        echo "reference_data_dir=$REFDATA" >> ${config_file}
        echo "reference_fasta=${ref_fasta}" >> ${config_file}
        echo "repeat_mask_file=${repeat_mask_bed}" >> ${config_file}
        echo "cutadapt_config_file=${cutadapt_config}" >> ${config_file}
        echo "cutadapt=${cutadapt_path}" >> ${config_file}
        echo "blat=${blat_path}" >> ${config_file}
        echo "jellyfish=${jellyfish_path}" >> ${config_file}
        echo "kmer_size=${kmer_size}" >> ${config_file}


        # Run BreaKmer
        python $BREAKMER -s ${indel_size} ${"-g " + gene_list} ${config_file}

        # Keep BreaKmer log file but discard gfServer log since it's not informative
        mv "${output_dir}/log.txt" "${output_log}"
        rm ${output_dir}/output/gfserver*.log

        # BreaKmer creates $OUTPUT/output directory structure, and this only saves
        # 'output' subdirectory that contains results
        mv ${output_dir} temp && mv temp/output ${output_dir}

        # Concatenate SV results
        if [[ `ls -1 ${output_dir}/*.out | grep "svs.out" | wc -l` -gt 0 ]]; then
            cat ${output_dir}/*_svs.out | grep ^gene | sort -u > ${sample_name}_svs.out
            cat ${output_dir}/*_svs.out | grep -v ^gene | sort -k1,1 >> ${sample_name}_svs.out
        else
            touch ${sample_name}_svs.out
        fi
    }

    output {
        File output_summary = "${output_dir}/${sample_name}_summary.out"
        File output_sv = "${sample_name}_svs.out"
        File breakmer_config = "${config_file}"
        File breakmer_log = "${output_log}"
    }

    runtime {
        docker: docker
        memory: select_first([memory_gb, 20]) + " GB"
        cpu: select_first([cpu, 1])
        disks: "local-disk " + disk_size + " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        maxRetries: select_first([max_retries, 1])
        bootDiskSizeGb: 10
    }
}

task GatherBreakmerResults {
    String sample_name
    Array[File] summary_files
    Array[File] sv_files

    command <<<

        # cosolidate summary files
        cat ${sep=" " summary_files} | grep ^Target | sort -u | cut -f1-6 > ${sample_name}_summary.out
        cat ${sep=" " summary_files} | grep -v ^Target | \
        awk 'BEGIN{OFS="\t"}
             {contig[$1]+=$2; tot_val[$1]+=$3; indel[$1]+=$4; sv[$1]+=$5; trl[$1]+=$6}
             END{for(i in contig){print i,contig[i],tot_val[i],indel[i],sv[i],trl[i]}}' |  \
        sort -k1,1 >> ${sample_name}_summary.out

        # concatenate SV results
        cat ${sep=" " sv_files} | grep ^gene | sort -u > ${sample_name}_sv.out
        cat ${sep=" " sv_files} | grep -v ^gene | sort -k1,1 >> ${sample_name}_sv.out

        # total SV detected
        cat ${sep=" " sv_files} | grep -v ^gene | wc -l > total_sv_num.txt
    >>>
    output {
        File output_summary = "${sample_name}_summary.out"
        File output_sv = "${sample_name}_sv.out"
        Int total_sv_detected = read_int("total_sv_num.txt")
    }
    runtime {
        docker: "ubuntu:16.04"
        memory: "2 GB"
        disks: "local-disk 10 HDD"
        preemptible: 3
        maxRetries: 1
    }
}
