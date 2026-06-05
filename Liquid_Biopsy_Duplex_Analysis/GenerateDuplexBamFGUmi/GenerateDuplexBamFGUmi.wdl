version 1.0
import "https://raw.githubusercontent.com/broadinstitute/TAG-public/f4554f714295ff2962136860195e7fb5727aa959/checkBaitSetName/checkBaitSetName.wdl" as checkBaitSetName

workflow GenerateDuplexBamFGUmi {
	input {
		String basename
		File input_bam
		File input_bam_idx
		File reference_fasta
		File reference_dict
		File reference_fasta_idx
		File reference_pac
		File reference_amb
		File reference_ann
		File reference_bwt
		File reference_sa
		File reference_0123
		File target_intervals
		File bait_intervals
		File known_indels
        File known_indels_index
        File variant_eval_gold_standard
        File variant_eval_gold_standard_index
        File dbsnp
        File dbsnp_index

		Boolean copy_umi_from_read_name = true
        Boolean run_bqsr = false
		String fgumi_docker = "us.gcr.io/tag-team-160914/fgumi:0.2.0"

		# parameters
		String bait_set
   		Boolean fail_on_intervals_mismatch = true
   		Int minimum_base_quality
   		Int allowable_umi_distance
   		String minimum_consensus_reads
   		String min_reads
   		Float frac_Ns
   		String? consensus_extra_filter_args
   		Int num_clip_bases_five_prime
   		Int? num_clip_bases_three_prime
		File process_duplex_coverage_rscript
		Int bwa_threads
		String bwa_path
		String bwa_docker
	}

	call checkBaitSetName.compareBaitSetName as checkBaitSetName {
      input:
         bait_set = bait_set,
         bait_intervals = bait_intervals,
         target_intervals = target_intervals,
         fail_task = fail_on_intervals_mismatch
   }
   if(copy_umi_from_read_name){
      call CopyUmiTask {
         input: 
         bam_file = input_bam,
         bam_index = input_bam_idx,
         base_name = basename
      }
   }
   call FGUmiSort {
	  input:
		basename = basename,
		input_bam = select_first([CopyUmiTask.umi_extracted_bam, input_bam]),
		input_bam_idx = select_first([CopyUmiTask.umi_extracted_bam_index, input_bam_idx]),
		docker_override = fgumi_docker
   }

   call FGUmiTask as FGUmi {
	  input:
		 basename = basename,
		 input_bam = FGUmiSort.output_bam,
		 target_intervals = target_intervals,
		 minimum_base_quality = minimum_base_quality,
		 allowable_umi_distance = allowable_umi_distance,
		 minimum_consensus_reads = minimum_consensus_reads,
		 min_reads = min_reads,
		 frac_Ns = frac_Ns,
		 consensus_extra_filter_args = consensus_extra_filter_args,
		 docker_override = fgumi_docker
   }

   call GetBwaVersion {
		input:
			docker=bwa_docker,
			bwa_path = bwa_path
	}

   call AlignDuplexConsensusReads {
		input:
		duplex_consensus_bam = FGUmi.duplex_output_bam,
		base_name = basename,
		reference_fasta = reference_fasta,
		reference_fasta_index = reference_fasta_idx,
		reference_dict = reference_dict,
		reference_pac = reference_pac,
		reference_amb = reference_amb,
		reference_ann = reference_ann,
		reference_bwt = reference_bwt,
		reference_sa = reference_sa,
		reference_0123 = reference_0123,
		threads = bwa_threads,
		docker_override = bwa_docker
	}

	call MergeBamAlignmentTask {
		input:
		reference = reference_fasta,
		reference_index = reference_fasta_idx,
		reference_pac = reference_pac,
		reference_amb = reference_amb,
		reference_ann = reference_ann,
		reference_bwt = reference_bwt,
		reference_sa = reference_sa,
		reference_dict = reference_dict,
		duplex_consensus_reads_bam = FGUmi.duplex_output_bam,
		aligned_duplex_reads_bam = AlignDuplexConsensusReads.aligned_bam,
		base_name = basename,
		bwa_path = bwa_path,
		bwa_version = GetBwaVersion.version,
		bwa_threads = bwa_threads
	}

	call ClipBam {
		input:
		bam_file = MergeBamAlignmentTask.output_bam,
		bam_index = MergeBamAlignmentTask.output_bam_index,
		base_name = basename,
		num_clip_bases_five_prime = num_clip_bases_five_prime,
		num_clip_bases_three_prime = num_clip_bases_three_prime,
		reference = reference_fasta,
		reference_index = reference_fasta_idx
	}

    if(run_bqsr){
        call BQSRWithoutBinning as BQSRDuplex {
          input:
           bam_file = ClipBam.output_bam,
           bam_index = ClipBam.output_bam_index,
           reference = reference_fasta,
           reference_index = reference_fasta_idx,
           reference_dict = reference_dict,
           dbsnp = dbsnp,
           dbsnp_index= dbsnp_index,
           known_indels = known_indels,
           known_indels_index = known_indels_index,
           variant_eval_gold_standard = variant_eval_gold_standard,
           variant_eval_gold_standard_index = variant_eval_gold_standard_index,
           base_name = basename
        }
    }

    call CollectDuplexDepthOfCoverage {
      input:
         interval_list = target_intervals,
         reference = reference_fasta,
         reference_index = reference_fasta_idx,
         reference_dict = reference_dict,
         bam_file = select_first([BQSRDuplex.output_bam, ClipBam.output_bam]),
         bam_index = select_first([BQSRDuplex.output_bam_index, ClipBam.output_bam_index]),
         base_name = basename,
         extra_arguments = "--countType COUNT_FRAGMENTS_REQUIRE_SAME_BASE -allowPotentiallyMisencodedQuals",
         process_duplex_coverage_rscript = process_duplex_coverage_rscript
   }

   call calculateDuplexDepthMetrics {
	  input:
		 depth_txt = CollectDuplexDepthOfCoverage.depth_txt,
		 duplex_depth_of_coverage = CollectDuplexDepthOfCoverage.duplex_depth_of_coverage
   }


   output {
		File duplex_output_bam = select_first([BQSRDuplex.output_bam, ClipBam.output_bam])
        File duplex_output_bam_index = select_first([BQSRDuplex.output_bam_index, ClipBam.output_bam_index])
		File duplex_family_sizes = FGUmi.duplex_family_sizes
		File duplex_yield_metrics = FGUmi.duplex_yield_metrics
		File family_sizes = FGUmi.family_sizes
		File umi_counts = FGUmi.umi_counts

        Int mean_duplex_depth = calculateDuplexDepthMetrics.mean_duplex_depth
        Float duplex_depth_above_500x = calculateDuplexDepthMetrics.duplex_depth_above_500x
        Float duplex_depth_above_1000x = calculateDuplexDepthMetrics.duplex_depth_above_1000x
        Float duplex_depth_above_1500x = calculateDuplexDepthMetrics.duplex_depth_above_1500x
        Float duplex_depth_above_2000x = calculateDuplexDepthMetrics.duplex_depth_above_2000x
   }
}

task FGUmiTask {
	input {
	String basename
	File input_bam
	File target_intervals
	Int minimum_base_quality
	Int allowable_umi_distance
	String minimum_consensus_reads
	String min_reads
	Float frac_Ns
	String? consensus_extra_filter_args

	String? docker_override
	String docker = select_first([docker_override, "us.gcr.io/tag-team-160914/fgumi:0.2.0"])
	Int preemptible_count = 2
	Int maxRetries = 1
	Int threads = 4
	Float mem = 32.0
	Int diskgb_buffer = 100

	}
	Int diskGb = ceil(size(input_bam,"GB")*5) + diskgb_buffer
	Int command_mem = ceil((mem * 1000 - 100) / threads)

	command {
		set -e

		mkdir -p /tmp

		fgumi group --threads ~{threads} --queue-memory ~{command_mem} \
		--input ~{input_bam} \
		--output ~{basename}.umi_grouped.bam \
		--strategy paired 

		fgumi duplex --threads ~{threads} --queue-memory ~{command_mem} \
		--input ~{basename}.umi_grouped.bam \
		--output ~{basename}.unfiltered_duplex.bam \
		--min-reads ~{min_reads}

		fgumi filter --threads ~{threads} --queue-memory ~{command_mem} \
		--input ~{basename}.unfiltered_duplex.bam \
		--output ~{basename}.duplex.bam \
		--min-reads ~{minimum_consensus_reads} \
		--min-base-quality ~{minimum_base_quality} \
		--max-no-call-fraction ~{frac_Ns} \
		~{consensus_extra_filter_args}

		fgumi duplex-metrics \
		--input ~{basename}.umi_grouped.bam \
		--intervals ~{target_intervals} \
		--output ~{basename}

		ls ./**
	}

	runtime {
		docker: docker
		disks: "local-disk " + diskGb + " HDD"
		memory: mem + " GB"
		cpu: threads
		maxRetries: maxRetries
		preemptible: preemptible_count
	}

	output {
		File duplex_output_bam = "~{basename}.duplex.bam"
		File duplex_family_sizes = "~{basename}.duplex_family_sizes.txt"
    	File duplex_yield_metrics = "~{basename}.duplex_yield_metrics.txt"
    	File family_sizes = "~{basename}.family_sizes.txt"
    	File umi_counts = "~{basename}.umi_counts.txt"
	}

}

task FGUmiSort {
	input {
	String basename
	File input_bam
	File input_bam_idx
	Boolean use_samtools = false

	String? docker_override
	String docker = select_first([docker_override, "us.gcr.io/tag-team-160914/fgumi:0.2.0"])
	Int preemptible_count = 2
	Int maxRetries = 1
	Int threads = 4
	Float mem = 64.0
	Int diskgb_buffer = 100
	String? sort_options
	}
	Int diskGb = ceil(size(input_bam,"GB")*5) + diskgb_buffer
	Int command_mem = ceil((mem * 1000 - 100) / threads / 2)

	command {
		set -e

		if [[ "~{use_samtools}" == "true" ]]; then
			samtools sort -@ ~{threads} -m ~{command_mem}M \
			--template-coordinate \
			-o ~{basename}.coordinate_sorted.bam ~{sort_options}\
			~{input_bam}
		else
			fgumi sort --threads ~{threads} --max-memory ~{command_mem}M \
			--input ~{input_bam} \
			--order template-coordinate \
			--output ~{basename}.coordinate_sorted.bam \
			~{sort_options}
		fi

	}

	runtime {
		docker: docker
		disks: "local-disk " + diskGb + " HDD"
		memory: mem + " GB"
		cpu: threads
		maxRetries: maxRetries
		preemptible: preemptible_count
	}

	output {
		File output_bam = "~{basename}.coordinate_sorted.bam"
	}

}

task CopyUmiTask {
	input {
    String bloodbiopsydocker = "us.gcr.io/tag-team-160914/liquidbiopsy:0.0.4.5"
    String base_name
    String? fgbio_override
    File bam_file
    File bam_index
    Boolean? remove_umi_from_read_name = true
    
    Int preemptible = 2
    Int maxRetries = 1
    Int? disk_pad
    Float? extra_mem
    Int cpu = 4
	}
    Int disk_size = ceil(size(bam_file, "GB") * 5) + select_first([disk_pad,0])
    Int compute_mem = ceil(mem) * 1000 - 500
    Float mem = 25 + select_first([extra_mem, 0])

    command {
        export FGBIO_LOCAL_JAR=~{default="/usr/fgbio-2.0.2.jar" fgbio_override}

        ln -vs ~{bam_file} ~{base_name}_input.bam
        ln -vs ~{bam_index} ~{base_name}_input.bai

        java -Xmx~{compute_mem}m -jar $FGBIO_LOCAL_JAR \
        CopyUmiFromReadName \
        -i ~{base_name}_input.bam \
        -o ~{base_name}.bam \
        --remove-umi ~{remove_umi_from_read_name}
    }
    
    output {
        File umi_extracted_bam = "~{base_name}.bam"
        File umi_extracted_bam_index = "~{base_name}.bai"
    }
    
    runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      maxRetries: maxRetries
      preemptible: preemptible
      cpu: cpu
  	}

}

task AlignDuplexConsensusReads {
	input {
		File duplex_consensus_bam
		String base_name
		File reference_fasta
		File reference_fasta_index
		File reference_dict
		File reference_pac
		File reference_amb
		File reference_ann
		File reference_bwt
		File reference_sa
		File reference_0123

		Int memory = 64
		Int disk_size = 200
		Int threads = 4
		Int cpu_cores = 1
		String? docker_override
		Int preemptible_attempts = 2
	}
		String output_bam_name = base_name + ".duplex_consensus.aligned_tmp.bam"
		String ref_name = basename(reference_fasta)

	command {
		set -e

		ln -vs ~{reference_fasta} ~{ref_name}
		ln -vs ~{reference_fasta_index} ~{ref_name}.fai
		ln -vs ~{reference_dict} ~{ref_name}.dict
		ln -vs ~{reference_bwt} ~{ref_name}.bwt.2bit.64
		ln -vs ~{reference_pac} ~{ref_name}.pac
		ln -vs ~{reference_amb} ~{ref_name}.amb
		ln -vs ~{reference_ann} ~{ref_name}.ann
		ln -vs ~{reference_sa} ~{ref_name}.sa
		ln -vs ~{reference_0123} ~{ref_name}.0123

		samtools fastq ~{duplex_consensus_bam} | \
		/opt/bwa-mem2/bwa-mem2 mem -K 100000000 -t ~{threads} -p -Y ~{ref_name} - > ~{base_name}.duplex_consensus.aligned_tmp.sam
		
		# This is to avoid the error of piping happening in multi-threading
		samtools view ~{base_name}.duplex_consensus.aligned_tmp.sam -S -b -o ~{output_bam_name}
	}

	output {
		File aligned_bam = "~{base_name}.duplex_consensus.aligned_tmp.bam"
	}

	runtime {
		memory: memory + " GB"
		docker: select_first([docker_override, "us.gcr.io/tag-public/bwa_mem2:latest"]) 
		disks: "local-disk " + disk_size + " HDD"
		cpu: cpu_cores
		preemptible: preemptible_attempts
	}
}

task MergeBamAlignmentTask {
	input {
	String bloodbiopsydocker = "us.gcr.io/tag-team-160914/liquidbiopsy:0.0.4.6"
	File? picard_override
	File reference
	File reference_index
	File reference_pac
	File reference_amb
	File reference_ann
	File reference_bwt
	File reference_sa
	File reference_dict
	File duplex_consensus_reads_bam
	File aligned_duplex_reads_bam
	String base_name
	String bwa_path
	String bwa_version
	Int bwa_threads
	Int? preemptible_attempts
	Int? memory
	Int disk_pad = 100
	}
	Int ref_size = ceil(size(reference, "GB") + size(reference_dict, "GB"))
	Int disk_size = ceil(size(duplex_consensus_reads_bam, "GB") * 15) + ref_size + disk_pad
	Int mem = select_first([memory, 32])

	command {
		set -e

		export PICARD_LOCAL_JAR=~{default="/usr/picard.jar" picard_override}

		java -jar -Xmx6144m $PICARD_LOCAL_JAR SortSam \
			INPUT=~{duplex_consensus_reads_bam} \
			OUTPUT=~{base_name}.sorted.unmapped.bam \
			SORT_ORDER=queryname

		java -jar -Xmx6144m $PICARD_LOCAL_JAR SortSam \
			INPUT=~{aligned_duplex_reads_bam} \
			OUTPUT=~{base_name}.sorted.mapped.bam \
			SORT_ORDER=coordinate

		java -jar -Xmx6144m $PICARD_LOCAL_JAR MergeBamAlignment \
			REFERENCE_SEQUENCE=~{reference} \
			UNMAPPED_BAM=~{base_name}.sorted.unmapped.bam \
			ALIGNED_BAM=~{base_name}.sorted.mapped.bam \
			OUTPUT=~{base_name}.merged.bam \
			PROGRAM_RECORD_ID="bwamem2" \
			PROGRAM_GROUP_VERSION="~{bwa_version}" \
			PROGRAM_GROUP_COMMAND_LINE="~{bwa_path} mem -K 100000000 -t ~{bwa_threads}" \
			PROGRAM_GROUP_NAME="bwamem2" \
		SORT_ORDER=coordinate \
		CREATE_INDEX=true \
			TMP_DIR=.

	}
	runtime {
		docker: bloodbiopsydocker
		disks: "local-disk " + disk_size + " HDD"
		memory: mem + "GB"
		cpu: "16"
		maxRetries: 3
		preemptible: select_first([preemptible_attempts, 10])
	}
	output {
		File output_bam = "~{base_name}.merged.bam"
		File output_bam_index = "~{base_name}.merged.bai"
	}
}

task GetBwaVersion {
	input {
		String docker
		String bwa_path
		Int? preemptible_attempts
	}
	command {
		~{bwa_path} version | tail -n1
	}
	runtime {
		docker: docker
		memory: "1 GB"
		maxRetries: 3
		preemptible: select_first([preemptible_attempts, 10])
	}
	output {
		String version = read_string(stdout())
	}
}

#want to clip the end of the fragments
#because they are lower quality
#also want to filter for MQ < 60
task ClipBam {
	input {
   String bloodbiopsydocker = "us.gcr.io/broad-dsde-methods/liquidbiopsy:0.0.3.5"
   File? picard_override
   File? fgbio_override
   File bam_file
   File bam_index
   String base_name
   Int num_clip_bases_five_prime
   File reference
   File reference_index
   Int? num_clip_bases_three_prime
   Int? preemptible_attempts
   Int disk_pad = 100
   Int? memory
   Int min_mapq = 60
	}
   Int ref_size = ceil(size(reference, "GB") +  size(reference_index, "GB"))
   Int disk_size = ceil(size(bam_file, "GB") * 5) + ceil(size(bam_index, "GB")) + ref_size + disk_pad
   Int mem = select_first([memory, 10])
   Int compute_mem = mem * 1000 - 1000

   command {

      set -e

      export FGBIO_LOCAL_JAR=~{default="/usr/fgbio-1.0.0.jar" fgbio_override}
      export PICARD_LOCAL_JAR=~{default="/usr/picard.jar" picard_override}

      java -Xmx~{compute_mem}m -jar $FGBIO_LOCAL_JAR ClipBam \
         -i ~{bam_file} \
         -o ~{base_name}.clipped.bam \
         -c "Hard" \
         --ref ~{reference} \
         --read-one-five-prime ~{num_clip_bases_five_prime} \
         --read-two-five-prime ~{num_clip_bases_five_prime} \
         ~{"--read-one-three-prime " + num_clip_bases_three_prime} \
         ~{"--read-two-three-prime " + num_clip_bases_three_prime}

      samtools view -hb -q ~{min_mapq} ~{base_name}.clipped.bam -o ~{base_name}.filtered.bam

      java -jar $PICARD_LOCAL_JAR BuildBamIndex \
         INPUT=~{base_name}.filtered.bam  \
         OUTPUT=~{base_name}.filtered.bai

   }
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 10])
   }
   output {
      File output_bam = "~{base_name}.filtered.bam"
      File output_bam_index = "~{base_name}.filtered.bai"
   }
}

task BQSRWithoutBinning {
    input {
   String bloodbiopsydocker = "us.gcr.io/broad-dsde-methods/liquidbiopsy:0.0.3.5"
   File? gatk_override
   File bam_file
   File bam_index
   File reference
   File reference_index
   File reference_dict
   File dbsnp
   File dbsnp_index
   File known_indels
   File known_indels_index
   File variant_eval_gold_standard
   File variant_eval_gold_standard_index
   String base_name
   Int? preemptible_attempts
   Int? memory
   Int disk_pad = 50
    }
   Int ref_size = ceil(size(reference, "GB") + size(reference_index, "GB") + size(reference_dict, "GB"))
   Int disk_size = ceil(size(bam_file, "GB") * 5) + ceil(size(bam_index, "GB")) + ref_size + disk_pad
   Int mem = select_first([memory, 5])
   Int compute_mem = mem * 1000 - 500

   command {
      set -e

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      gatk BaseRecalibrator \
         -R ~{reference} \
         -I ~{bam_file} \
         --known-sites ~{dbsnp} \
         --known-sites ~{known_indels} \
         --known-sites ~{variant_eval_gold_standard} \
         -O ~{base_name}.recalibration_report.grp

      gatk ApplyBQSR \
         -R ~{reference} \
         -I ~{bam_file} \
         -bqsr ~{base_name}.recalibration_report.grp \
         -O ~{base_name}.bqsr.bam
   }
   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 2])
   }
   output {
      File output_bam = "~{base_name}.bqsr.bam"
      File output_bam_index = "~{base_name}.bqsr.bai"
      File output_recal_report_grp = "~{base_name}.recalibration_report.grp"
   }
}

task CollectDuplexDepthOfCoverage {
   input {
    File interval_list
    File reference
    File reference_index
    File reference_dict
    File bam_file
    File bam_index
    File process_duplex_coverage_rscript
    String base_name
    String? gatk3_override
    String? extra_arguments
    Int? preemptible_attempts
    Int? memory
    String docker = "us.gcr.io/tag-public/tag-tools:1.0.0"
    Int disk_pad = 50
   }
   Int ref_size = ceil(size(reference, "GB") + size(reference_index, "GB") + size(reference_dict, "GB"))
   Int disk_size = ceil(size(bam_file, "GB") * 2) + ceil(size(bam_index, "GB")) + ref_size + disk_pad
   Int mem = select_first([memory, 15])
   Int compute_mem = mem * 1000 - 1000
   
   command {
      set -e

      export GATK_LOCAL_JAR=~{default="/usr/tag/GATK36.jar" gatk3_override}

      # Calculate tumor depth over the panel
      # only count fragments with same base.
      java -Xmx~{compute_mem}m -jar $GATK_LOCAL_JAR -T DepthOfCoverage \
         -L ~{interval_list} \
         -I ~{bam_file} \
         -R ~{reference} \
         -o ~{base_name}.depth \
         --omitPerSampleStats \
         --printBaseCounts \
         ~{extra_arguments}

      Rscript -e 'source("~{process_duplex_coverage_rscript}"); generateDepthFigures("~{base_name}", "~{base_name}.depth")'
   }

   runtime {
      docker: docker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + "GB"
      maxRetries: 3
      preemptible: select_first([preemptible_attempts, 2])
   }

   output {
      File depth_txt = "~{base_name}.depth.txt"
	  File duplex_depth_of_coverage = "~{base_name}.depth"
   }
}

task calculateDuplexDepthMetrics {
    input {
        File depth_txt
        File duplex_depth_of_coverage
    }

    command {
        set -e
        python <<CODE

      import pandas as pd

      df = pd.read_csv("~{depth_txt}", delim_whitespace=True)

      def writeFile(value, filename):
         f = open(filename + ".txt", 'w')
         f.write(str(int(float(value))))
         f.close()

      for col in df.columns:
         writeFile(df[col].loc[0], col)

      def writeDepthStatistic(depthValues, depthCutoff):
         f = open(str(depthCutoff) + ".txt", 'w')
         f.write(str(depthValues[depthValues.Total_Depth >= depthCutoff].count().Total_Depth / depthValues.count().Total_Depth))
         f.close()

      depthOfCoverageByLocus = pd.read_csv("~{duplex_depth_of_coverage}", delimiter = '\t')

      writeDepthStatistic(depthOfCoverageByLocus, 500)
      writeDepthStatistic(depthOfCoverageByLocus, 1000)
      writeDepthStatistic(depthOfCoverageByLocus, 1500)
      writeDepthStatistic(depthOfCoverageByLocus, 2000)

      CODE
	}

	output {
      Int mean_duplex_depth = read_int('meanDuplexDepth.txt')
      Float duplex_depth_above_500x = read_float('500.txt')
      Float duplex_depth_above_1000x = read_float('1000.txt')
      Float duplex_depth_above_1500x = read_float('1500.txt')
      Float duplex_depth_above_2000x = read_float('2000.txt')
	}

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/liquidbiopsy:0.0.3.7"
        memory: "2 GB"
        cpu: "1"
        disks: "local-disk 25 HDD"
    }
}