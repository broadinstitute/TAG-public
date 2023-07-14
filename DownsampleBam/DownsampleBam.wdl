workflow DownsampleBam{
    String output_name
    File? target_intervals
    File? bait_intervals
    Boolean is_wgs 
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fai
    
    if(is_wgs){
        call CollectWgsMetrics as CollectWgsMetricsBeforeDownsample{
  			input: bam = input_bam,
        	   	   bam_index = input_bam_index,
        	   	   output_name = output_name + ".pre_downsample",
                   ref_fasta = ref_fasta,
                   interval_list = target_intervals,
                   ref_fai = ref_fai
    	}
        call GetMeanCoverage {
        	input: wgs_metrics = CollectWgsMetricsBeforeDownsample.metrics
        }
    }
    if(!is_wgs){
        call CollectHsMetrics as CollectHsMetricsBeforeDownsample{
  			input: bam = input_bam,
        	   	   bam_index = input_bam_index,
                   target_intervals = target_intervals,
                   bait_intervals = bait_intervals,
        	   	   output_name = output_name + ".pre_downsample"
    	}
        call GetMeanTargetCoverage {
    		input: hs_metrics = CollectHsMetricsBeforeDownsample.metrics
    	}
    }
    call MarkDuplicates {
    	input: bam = input_bam,
        	   bam_index = input_bam_index
    }
    call DownsampleSam {
    	input: output_name=output_name,
        	   bam = MarkDuplicates.output_bam,
               bam_index = MarkDuplicates.output_bam_index,
               coverage = select_first([GetMeanTargetCoverage.target_coverage, GetMeanCoverage.coverage])
    }
    if(is_wgs){
    	call CollectWgsMetrics {
  			input: bam = DownsampleSam.output_bam,
        	   	   bam_index = DownsampleSam.output_bam_index,
        	   	   output_name = output_name,
                   interval_list = target_intervals,
                   ref_fasta = ref_fasta,
                   ref_fai = ref_fai
    	}
        call GetMeanCoverage as GetDSCoverage {
    		input: wgs_metrics = CollectWgsMetrics.metrics
    	}
    }
    if(!is_wgs){
    	call CollectHsMetrics {
  			input: bam = DownsampleSam.output_bam,
        	   	   bam_index = DownsampleSam.output_bam_index,
                   target_intervals = target_intervals,
                   bait_intervals = bait_intervals,
        	   	   output_name = output_name
    	}
        call GetMeanTargetCoverage as GetDSTargetCoverage{
    		input: hs_metrics = CollectHsMetrics.metrics
    	}
    }
    
    output {
    	File metrics = select_first([CollectWgsMetrics.metrics, CollectHsMetrics.metrics])
        File pre_downsample_metrics = select_first([CollectWgsMetricsBeforeDownsample.metrics, CollectHsMetricsBeforeDownsample.metrics])
    	File ds_bam = DownsampleSam.output_bam
        File ds_bam_index = DownsampleSam.output_bam_index
        File md_metrics = MarkDuplicates.duplicate_metrics
        Float? pre_downsampled_mean_coverage = GetMeanCoverage.coverage
        Float? pre_downsampled_mean_target_coverage = GetMeanTargetCoverage.target_coverage
        Float? downsampled_mean_coverage = GetDSCoverage.coverage
        Float? downsampled_mean_target_coverage = GetDSTargetCoverage.target_coverage
    }

}

task GetMeanTargetCoverage {

    File hs_metrics
    Int? extra_disk
    Int disk_size = 50 + select_first([extra_disk,0])

	command <<<
    	coverage="$(tail -n +8 ${hs_metrics} | head -n1 | cut -f34)"
		echo $coverage
	>>>

	runtime {
    	preemptible: 3
    	docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:7bc64948a0a9f50ea55edb8b30c710943e44bd861c46a229feaf121d345e68ed"
	 	memory: "2 GB"
     	cpu: "1"
		disks: "local-disk " + disk_size+ " HDD"
	}
	output {
		 Float target_coverage = read_float(stdout())
	}
}

task GetMeanCoverage {

    File wgs_metrics
    Int? extra_disk
    Int disk_size = 50 + select_first([extra_disk,0])

	command <<<
     index="$(grep MEAN_COVERAGE ${wgs_metrics} | awk -v name='MEAN_COVERAGE' '{for (i=1;i<=NF;i++) if ($i==name) print i; exit}' )"
		 coverage="$(grep -v '#' ${wgs_metrics} | awk -v i=$index 'FNR == 3 {print $i}')"
		 echo $coverage
	>>>

	runtime {
    	preemptible: 3
    	docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:7bc64948a0a9f50ea55edb8b30c710943e44bd861c46a229feaf121d345e68ed"
	 	memory: "2 GB"
     	cpu: "1"
		disks: "local-disk " + disk_size+ " HDD"
	}
	output {
		 Float coverage = read_float(stdout())
	}
}

task CollectWgsMetrics {
    File bam
    File bam_index 
    File ref_fasta
    File ref_fai
    String output_name
    File? interval_list
    Int? extra_disk
    Float? extra_mem
    Int disk_size = ceil(size(bam,"GB") * 1.5) + 20 + select_first([extra_disk,0])
    Float memory = 5 + select_first([extra_mem, 0])

	command {
		/gatk/gatk CollectWgsMetrics -I ${bam} -O ${output_name}.wgs_metrics -R ${ref_fasta} ${'--INTERVALS ' + interval_list}
	}

	runtime {
         docker: "broadinstitute/gatk:4.2.5.0"
         memory: memory + " GB"
         cpu: "1"
         disks: "local-disk " + disk_size + " HDD"
         preemptible: 2
    }

	output {
		File metrics = "${output_name}.wgs_metrics"
	}
}

task CollectHsMetrics {
    File bam
    File bam_index
    String output_name
    File bait_intervals
    File target_intervals
    Int? extra_disk
    Float? extra_mem
    Int disk_size = ceil(size(bam,"GB") * 1.5) + 20 + select_first([extra_disk,0])
    Float memory = 5 + select_first([extra_mem, 0])

	command {
		/gatk/gatk CollectHsMetrics -I ${bam} -O ${output_name}.hs_metrics -BI ${bait_intervals} -TI ${target_intervals}
	}

	runtime {
         docker: "broadinstitute/gatk:4.2.5.0"
         memory: memory + " GB"
         cpu: "1"
         disks: "local-disk " + disk_size + " HDD"
         preemptible: 2
    }

	output {
		File metrics = "${output_name}.hs_metrics"
	}
}

task MarkDuplicates {
   File bam
   File bam_index
   File? ref_fasta
   File? ref_fasta_index
   File? ref_fasta_dict
   Int? extra_disk
   Int disk_size = ceil(size(bam,"GB") * 1.5) + 20 + select_first([extra_disk,0])
   Int? preemptible_attempts
   Int? memory
   Int mem = select_first([memory, 8])
   String out_basename = basename(bam, ".bam")

   command {
   
      /gatk/gatk MarkDuplicates \
         INPUT=${bam} \
         ${if defined("ref_fasta") then "REFERENCE_SEQUENCE=" + ref_fasta else ""} \
         REMOVE_DUPLICATES=true \
         OUTPUT=${out_basename}.bam \
         M=${out_basename}.duplicate_metrics.txt \
         CREATE_INDEX=true
   }
   runtime {
      docker: "broadinstitute/gatk:4.2.5.0"
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      preemptible: select_first([preemptible_attempts, 3])
   }
   output {
      File output_bam = "${out_basename}.bam"
      File output_bam_index = "${out_basename}.bai"
      File duplicate_metrics = "${out_basename}.duplicate_metrics.txt"
   }
}
task DownsampleSam {
    File bam
    File bam_index
    File? ref_fasta
    File? ref_fasta_index
    File? ref_fasta_dict
    Float coverage
    Float? desired_coverage
    Float? desired_frac
    Int? random_seed
    String output_name
    Int? extra_disk
    Float? extra_mem 
    Float? downsampleFraction = if defined(desired_frac) then desired_frac else desired_coverage/coverage
    Int disk_size = ceil(size(bam,"GB") * 2.5 + 50) + select_first([extra_disk,0])
    Float memory = 7.5 + select_first([extra_mem, 0])

	command <<<
    	dsFrac="$(awk -v frac=${downsampleFraction} 'BEGIN{if(frac>1.0){print 1.0}else{print frac}}')"
		/gatk/gatk DownsampleSam --INPUT ${bam} --OUTPUT ${output_name}.ds.bam ${'--REFERENCE_SEQUENCE ' + ref_fasta} -P $dsFrac --CREATE_INDEX true ${"-R " + random_seed}
	>>>


	runtime {
         docker: "broadinstitute/gatk:4.2.5.0"
         memory: memory + " GB"
         cpu: "1"
         disks: "local-disk " + disk_size + " HDD"
    }

	output {
		File output_bam = "${output_name}.ds.bam"
		File output_bam_index = "${output_name}.ds.bai"
	}
}


