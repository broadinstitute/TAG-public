version 1.0

import "../../checkBaitSetName/checkBaitSetName.wdl" as checkBaitSetName
import "./mutect2.wdl" as m2


workflow MakeCallsFromConsensus {

   input {
      # dockers and override jars
      String bloodbiopsydocker 
      String reference_version
      Int preemptible

      # Mutect2 specific arguments
      String basename
      File tumor_bam
      File tumor_bam_idx
      String tumor_sample_name
      File? normal_bam
      File? normal_bam_idx
      String? normal_sample_name
      File? raw_tumor_bam
      File? raw_tumor_bam_idx
      File? raw_normal_bam
      File? raw_normal_bam_idx

      File target_intervals
      Boolean fail_on_intervals_mismatch
      File reference
      File reference_idx
      File reference_dict
      File? pon
      File? pon_idx
      File? pon_for_benchmarking
      File? pon_for_benchmarking_idx
      File? variants_for_contamination
      File? variants_for_contamination_idx
      File? gnomad
      File? gnomad_idx
      File? gnomad_for_benchmarking
      File? gnomad_for_benchmarking_idx
      Boolean? run_orientation_bias_mixture_model_filter
      Boolean run_ob_filter = select_first([run_orientation_bias_mixture_model_filter, false])
      Boolean FilterVariants
      Int scatter_count
      String? m2_extra_args
      String? m2_extra_filtering_args
      String m2_extra_filtering_args_or_default = select_first([m2_extra_filtering_args, ""])
      Boolean? compress_vcfs
      String mapping_filter_python_script
      File blastdb_nhr
      File blastdb_nin
      File blastdb_nsq
      String blastn_path

      Boolean? is_benchmark

      File? gatk_override
   
      # Fingerprinting haplotype database
      File haplotype_map
      Boolean? fingerprint_tumor_normal

      ## Use as a last resort to increase the disk given to every task in case of ill behaving data
      Int? emergency_extra_disk
      # This is added to every task as padding, should increase if systematically you need more disk for every call
      Int disk_pad = 10 + select_first([emergency_extra_disk,0])
      #if set to true, use the strand bias filter on the raw annotation
      #if set to false (default), uses the strand bias filter on the duplex annotations

      # Funcotator inputs
      Boolean filter_funcotations
      File? data_sources_tar_gz
      String? funcotator_extra_args
   }

   Int small_task_mem = 4
   Int small_task_disk = 100
   Int boot_disk_size = 12
   Int small_task_cpu = 2
   Int max_retries_or_default = 2

   Runtime standard_runtime = {"gatk_docker": bloodbiopsydocker, "gatk_override": gatk_override,
      "max_retries": max_retries_or_default, "preemptible": preemptible, "cpu": small_task_cpu,
      "machine_mem": small_task_mem * 1000, "command_mem": small_task_mem * 1000 - 500,
      "disk": small_task_disk + disk_pad, "boot_disk_size": boot_disk_size}
   
   if(select_first([fingerprint_tumor_normal, true])) {
      call CrosscheckFingerprints {
         input:
            bloodbiopsydocker = bloodbiopsydocker,
            preemptible = preemptible,
            disk_pad = disk_pad,
            input_bam = tumor_bam,
            second_input_bam = normal_bam,
            haplotype_map = haplotype_map
      }
      call FingerprintsResult{
         input:
            bloodbiopsydocker = bloodbiopsydocker,
            preemptible = preemptible,
            disk_pad = disk_pad,
            fingerprint_metrics = CrosscheckFingerprints.fingerprint_metrics,

      }
   }
   call checkBaitSetName.compareBaitSetName as checkBaitSetName {
      input:
         target_intervals = target_intervals,
         fail_task = fail_on_intervals_mismatch
   }
   # Collect Sequencing Artifact Metrics after deduplication by start and stop
   # positions (but not including UMIs).
   call CollectSequencingArtifactMetrics as ConsensusArtifactMetricsTumor {
      input:
         bam_file = tumor_bam,
         bam_idx = tumor_bam_idx,
         basename = tumor_sample_name,
         reference = reference,
         reference_idx = reference_idx,
         bloodbiopsydocker = bloodbiopsydocker,
         preemptible = preemptible,
         disk_pad = disk_pad
   }

   if (defined(normal_bam)) {
      call CollectSequencingArtifactMetrics as ConsensusArtifactMetricsNormal {
         input:
            bam_file = select_first([normal_bam]),
            bam_idx =  select_first([normal_bam_idx]),
            basename = normal_sample_name,
            reference = reference,
            reference_idx = reference_idx,
            bloodbiopsydocker = bloodbiopsydocker,
            preemptible = preemptible,
            disk_pad = disk_pad
      }
   }

   Boolean is_benchmark_or_default = select_first([is_benchmark, false])
   File? pon_to_use = if is_benchmark_or_default then pon_for_benchmarking else pon
   File? pon_to_use_idx = if is_benchmark_or_default then pon_for_benchmarking_idx else pon_idx

   File? gnomad_to_use = if is_benchmark_or_default then gnomad_for_benchmarking else gnomad
   File? gnomad_to_use_idx = if is_benchmark_or_default then gnomad_for_benchmarking_idx else gnomad_idx

   call m2.Mutect2 as M2Duplex {
      input:
          intervals = target_intervals,
          ref_fasta = reference,
          ref_fai = reference_idx,
          ref_dict = reference_dict,
          tumor_reads = tumor_bam,
          tumor_reads_index = tumor_bam_idx,
          normal_reads = normal_bam,
          normal_reads_index = normal_bam_idx,
          pon = pon_to_use,
          pon_idx = pon_to_use_idx,
          scatter_count = scatter_count,
          gnomad = gnomad_to_use,
          gnomad_idx = gnomad_to_use_idx,
          variants_for_contamination = variants_for_contamination,
          run_orientation_bias_mixture_model_filter = run_ob_filter,
          m2_extra_args = m2_extra_args,
          m2_extra_filtering_args = m2_extra_filtering_args_or_default,
          gatk_docker = bloodbiopsydocker,
          compress_vcfs = compress_vcfs
   }
   
   if(defined(M2Duplex.contamination_table)){
       call ExtractContam {
           input: bloodbiopsydocker = bloodbiopsydocker,
                  preemptible = preemptible,
                  contamination_table = select_first([M2Duplex.contamination_table])
       }
   }

   if (FilterVariants) {
   # Apply MRD mapping filter to duplex
   call RunMappingFilter {
      input:
         bloodbiopsydocker = bloodbiopsydocker,
         basename = basename,
         blastn_path = blastn_path,
         mapping_filter_python_script = mapping_filter_python_script,
         vcf_file = M2Duplex.filtered_vcf,
         vcf_idx = M2Duplex.filtered_vcf_idx,
         reference = reference,
         reference_idx = reference_idx,
         reference_dict = reference_dict, 
         blastdb_nhr = blastdb_nhr,
         blastdb_nin = blastdb_nin,
         blastdb_nsq = blastdb_nsq, 
         preemptible = preemptible,
         disk_pad = disk_pad
   }

   # Run additional filtering tasks
   call VariantFiltration {
      input: 
         bloodbiopsydocker = bloodbiopsydocker,
         basename = basename,
         mapping_filter_name = "mapping_filter",
         duplex_vcf_file = M2Duplex.filtered_vcf,
         duplex_vcf_idx = M2Duplex.filtered_vcf_idx,
         filter_vcf_file = RunMappingFilter.map_filtered_vcf,
         filter_vcf_idx = RunMappingFilter.map_filtered_vcf_idx,
         reference_fasta = reference,
         reference_fasta_idx = reference_idx,
         reference_dict = reference_dict, 
         preemptible = preemptible,
         filter_not_in_mask = true,
         disk_pad = disk_pad
   }
   }
   
   File filtered_vcf =  select_first([VariantFiltration.output_vcf, M2Duplex.filtered_vcf])
   File filtered_vcf_idx = select_first([VariantFiltration.output_vcf_idx, M2Duplex.filtered_vcf_idx])

   # Split VCFs by snps and indels
   call SplitVCFs {
      input:
         bloodbiopsydocker = bloodbiopsydocker,
         basename = basename,
         vcf_file =  M2Duplex.filtered_vcf,
         vcf_idx =  M2Duplex.filtered_vcf_idx,
         reference_fasta = reference,
         reference_fasta_idx = reference_idx,
         reference_dict = reference_dict, 
         preemptible = preemptible,
         disk_pad = disk_pad
   }

   # Create text file of variants
   call VariantsToTable {
      input: 
         bloodbiopsydocker = bloodbiopsydocker,
         vcf = SplitVCFs.output_snp_vcf, 
         vcf_idx = SplitVCFs.output_snp_vcf_idx,
         output_name = basename, 
         preemptible = preemptible,
         disk_pad = disk_pad
   } 

   # Create annotated maf file containing SNPs
   call m2.Funcotate as FuncotateMafSnps {
      input:
         ref_fasta = reference,
         input_vcf = SplitVCFs.output_snp_vcf,
         input_vcf_idx = SplitVCFs.output_snp_vcf_idx,
         case_id = tumor_sample_name,
         control_id = normal_sample_name,
         reference_version = reference_version,
         data_sources_tar_gz = data_sources_tar_gz,
         filter_funcotations = filter_funcotations,
         extra_args = funcotator_extra_args,
         use_gnomad = false, 
         output_format = "MAF",
         output_file_base_name = basename + "_snp", 
         compress = true,
         runtime_params = standard_runtime
   }

   # Create annotated maf file containing indels
   call m2.Funcotate as FuncotateMafIndels {
      input:
         ref_fasta = reference,
         input_vcf = SplitVCFs.output_indel_vcf,
         input_vcf_idx = SplitVCFs.output_indel_vcf_idx,
         case_id = tumor_sample_name,
         control_id = normal_sample_name,
         reference_version = reference_version,
         data_sources_tar_gz = data_sources_tar_gz,
         filter_funcotations = filter_funcotations,
         extra_args = funcotator_extra_args,
         use_gnomad = false, 
         output_format = "MAF",
         output_file_base_name = basename + "_indel", 
         compress = true,
         runtime_params = standard_runtime
   }

   Array[String] input_files = select_all([SplitVCFs.output_snp_vcf, SplitVCFs.output_indel_vcf, tumor_bam, normal_bam])

   # Create IGV session so that it is easy to examine results of pipeline
   call GenerateIGVSession {
      input: 
         bloodbiopsydocker = bloodbiopsydocker,
         preemptible = preemptible,
         input_files = input_files,
         file_name =  basename, 
         reference_version = reference_version
   }
   
   output {

      File igv_session = GenerateIGVSession.igv_session

      File filtered_vcf = VariantFiltration.output_vcf
      File filtered_vcf_idx = VariantFiltration.output_vcf_idx

      File filtered_snp_vcf = SplitVCFs.output_snp_vcf
      File filtered_snp_vcf_idx = SplitVCFs.output_snp_vcf_idx

      File filtered_indel_vcf = SplitVCFs.output_indel_vcf
      File filtered_indel_vcf_idx = SplitVCFs.output_indel_vcf_idx

      Int n_passing_snps = SplitVCFs.passing_SNP
      Int n_filtered_snps = SplitVCFs.filtered_SNP
      
      Int n_passing_mnps = SplitVCFs.passing_MNP
      Int n_filtered_mnps = SplitVCFs.filtered_MNP
      
      Int n_passing_indels = SplitVCFs.passing_INDEL
      Int n_filtered_indels = SplitVCFs.filtered_INDEL

      File? contamination_table = M2Duplex.contamination_table
      Float? contamination_fraction = ExtractContam.fracContam

      File funcotated_snp_maf = FuncotateMafSnps.funcotated_output_file
      File funcotated_indel_maf = FuncotateMafIndels.funcotated_output_file

      File variant_table = VariantsToTable.output_table
      File? fingerprint_metrics = CrosscheckFingerprints.fingerprint_metrics
      Int? expected_match = FingerprintsResult.expected_match

   }
}

# If pipeline is run in tumor/normal mode, ensure that they share the same fingerprint
task CrosscheckFingerprints {
   input {
      String bloodbiopsydocker
      Int preemptible
      Int disk_pad
      File input_bam
      File? second_input_bam
      File haplotype_map
      Int? memory

      Int disk_size = 20 + disk_pad
      Int mem = select_first([memory, 5])
   }
   command {
      set -e
      
      gatk CrosscheckFingerprints \
         --I ~{input_bam} \
         ~{"--I " + second_input_bam} \
         --EXIT_CODE_WHEN_MISMATCH 1 \
         --CROSSCHECK_BY SAMPLE \
         --EXPECT_ALL_GROUPS_TO_MATCH true \
         --HAPLOTYPE_MAP ~{haplotype_map} \
         --OUTPUT fingerprint_metrics.txt
    }
    runtime {
       docker: bloodbiopsydocker
       memory: mem + " GB"
       disks: "local-disk " + disk_size + " HDD"
       preemptible: preemptible
    }
    output {
       File fingerprint_metrics = "fingerprint_metrics.txt"
    }
}

task FingerprintsResult{
   input{
      File fingerprint_metrics
      String bloodbiopsydocker
      Int preemptible
      Int disk_pad
      Int disk_size = 20 + disk_pad
      Int? memory
      Int mem = select_first([memory, 5])
   }
   command{
      # Extract output from the fingerprint matrix

      python3 <<CODE 
      import pandas as pd

      fps_mtx = pd.read_csv("${fingerprint_metrics}",sep='\t',skiprows=6)
      # Get a brief output from crosscheck-result.txt
      # if all samples matched with each other 
      # it means they are very likely from the same individual
      with open('crosscheck-result.txt', 'w') as f:
         if 'EXPECTED_MATCH' in fps_mtx['RESULT'].unique():
               f.write('1')
         else:
               f.write('0')
      CODE
   }
   runtime {
       docker: bloodbiopsydocker
       memory: mem + " GB"
       disks: "local-disk " + disk_size + " HDD"
       preemptible: preemptible
    }
   output{
       Int expected_match = read_int('crosscheck-result.txt')
   }
}


# Collect Sequencing Artifact Metrics
# This is generally a good set of metrics to collect, but
# we will also need this later for the variant filtration
# step in FilterByOrientationBias
task CollectSequencingArtifactMetrics {
   input {
      String bloodbiopsydocker
      File bam_file
      File bam_idx
      File reference
      File reference_idx
      String? basename
      Int disk_pad
      Int preemptible
      Int ref_size = ceil(size(reference, "GB") + size(reference_idx, "GB"))
      Int? memory
      Int disk_size = 500
   }
   Int mem = select_first([memory, 5])
   Int compute_mem = mem * 1000 - 500

   command {
      set -e

      gatk --java-options "-Xmx4G" \
         CollectSequencingArtifactMetrics \
         --I ${bam_file} \
         --O "${basename}.gatk" \
         --R ${reference} \
         --VALIDATION_STRINGENCY LENIENT
   }
   runtime {
      docker: bloodbiopsydocker
      memory: mem + " GB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: preemptible
   }
   output {
      File pre_adapter_metrics = "${basename}.gatk.pre_adapter_detail_metrics"
   }
}


task RunMappingFilter {

   input {
      String bloodbiopsydocker
      String basename
      File vcf_file
      File vcf_idx
      File reference
      File reference_idx
      File reference_dict
      String mapping_filter_python_script
      String blastn_path
      File blastdb_nhr
      File blastdb_nin
      File blastdb_nsq
      Int preemptible
      Int disk_pad
      Int? memory
   }

   Int ref_size = ceil(size(reference, "GB") + size(reference_idx, "GB") + size(reference_dict, "GB"))
   Int blast_size = ceil(size(blastdb_nhr, "GB") + size(blastdb_nin, "GB") + size(blastdb_nsq, "GB"))
   Int mem = select_first([memory, 10])
   Int disk_size = ceil(size(vcf_file, "GB") * 2) + ref_size + blast_size + disk_pad

   command {

      set -e

      python2.7 ${mapping_filter_python_script} \
         --vcf ${vcf_file} \
         --outfile ${basename}.filtered.vcf \
         --reference_fasta ${reference} \
         --blastn ${blastn_path}

      bgzip "${basename}.filtered.vcf"
      tabix "${basename}.filtered.vcf.gz"

   }
   runtime {
      docker: bloodbiopsydocker
      preemptible: preemptible
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
   }
   output {
      File map_filtered_vcf = "${basename}.filtered.vcf.gz"
      File map_filtered_vcf_idx = "${basename}.filtered.vcf.gz.tbi"
   }
}


# This task filters out anything not included in the filter_vcf_file
task VariantFiltration {

   input {
      String bloodbiopsydocker
      String basename
      File? gatk_override
      String mapping_filter_name
      File duplex_vcf_file
      File duplex_vcf_idx
      File filter_vcf_file
      File filter_vcf_idx
      File reference_fasta
      File reference_fasta_idx
      File reference_dict
      Int preemptible
      Boolean filter_not_in_mask
      Int? memory
      Int disk_pad
   }

   Int mem = select_first([memory, 5])
   Int ref_size = ceil(size(reference_fasta, "GB") + size(reference_fasta_idx, "GB") + size(reference_dict, "GB"))
   Int disk_size = ceil(size(duplex_vcf_file, "GB") * 1.25) + ceil(size(filter_vcf_file, "GB") * 1.25) + ref_size + disk_pad

   command {
      set -e

      export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

      # Filter out variant with POPAF < 3.0 (this is phred scaled), which is
      # equivalent to filtering population allele frequency > 0.001.
      gatk VariantFiltration \
         -R ${reference_fasta} \
         -V ${duplex_vcf_file} \
         -O popaf.filtered.vcf.gz \
         --filter-expression "POPAF < 3.0" \
         --filter-name germline

      # Apply mapping filter
      gatk VariantFiltration \
         -R ${reference_fasta} \
         -V popaf.filtered.vcf.gz \
         -O ${basename}.${mapping_filter_name}.vcf.gz \
         --mask ${filter_vcf_file} \
         --filter-not-in-mask ${filter_not_in_mask} \
         --mask-name ${mapping_filter_name}

   }

   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      preemptible: preemptible
   }

   output {
      File output_vcf = "${basename}.${mapping_filter_name}.vcf.gz"
      File output_vcf_idx = "${basename}.${mapping_filter_name}.vcf.gz.tbi"
   }

}


task SplitVCFs {

   input {
      String bloodbiopsydocker
      String basename
      File? gatk_override
      File vcf_file
      File vcf_idx
      File reference_fasta
      File reference_fasta_idx
      File reference_dict
      Int preemptible
      Int? memory
      Int disk_pad
   }

   Int ref_size = ceil(size(reference_fasta, "GB") + size(reference_fasta_idx, "GB") + size(reference_dict, "GB"))
   Int disk_size = ceil(size(vcf_file, "GB") * 2) + ceil(size(vcf_idx)) + ref_size + disk_pad
   Int mem = select_first([memory, 16])

   command {

      set -e
      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      # Include all snp and mnp calls in the snp vcf
      gatk --java-options "-Xmx15G" SelectVariants \
         -V ${vcf_file} \
         -O "${basename}.snp.vcf.gz" \
         -R ${reference_fasta} \
         -select-type SNP \
         -select-type MNP

      #passing only SNPs 
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "${basename}.snp.vcf.gz" \
         -O "${basename}.snp.passing.vcf.gz" \
         -R ${reference_fasta} \
         --exclude-filtered \
         -select-type SNP

      passing_SNP="$(gatk CountVariants -V ${basename}.snp.passing.vcf.gz | tail -1)"
      echo "$passing_SNP" > passing_snp.txt
      
      #passing only MNPs
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "${basename}.snp.vcf.gz" \
         -O "${basename}.mnp.passing.vcf.gz" \
         -R ${reference_fasta} \
         --exclude-filtered \
         -select-type MNP
         
      passing_MNP="$(gatk CountVariants -V ${basename}.mnp.passing.vcf.gz | tail -1)"
      echo "$passing_MNP" > passing_mnp.txt

      #filtered only SNPs
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "${basename}.snp.vcf.gz" \
         -XL "${basename}.snp.passing.vcf.gz" \
         -O "${basename}.snp.filtered.vcf.gz" \
         -select-type SNP \
         -R ${reference_fasta} 

      filtered_SNP="$(gatk CountVariants -V ${basename}.snp.filtered.vcf.gz | tail -1)"
      echo "$filtered_SNP" > filtered_snp.txt
      
      #filtered only MNPs
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "${basename}.snp.vcf.gz" \
         -XL "${basename}.mnp.passing.vcf.gz" \
         -O "${basename}.mnp.filtered.vcf.gz" \
         -select-type MNP \
         -R ${reference_fasta} 

      filtered_MNP="$(gatk CountVariants -V ${basename}.mnp.filtered.vcf.gz | tail -1)"
      echo "$filtered_MNP" > filtered_mnp.txt

      gatk --java-options "-Xmx15G" SelectVariants \
         -V ${vcf_file} \
         -O "${basename}.indel.vcf.gz" \
         -R ${reference_fasta} \
         -select-type INDEL \
         -select-type MIXED

      #passing only 
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "${basename}.indel.vcf.gz" \
         -O "${basename}.indel.passing.vcf.gz" \
         -R ${reference_fasta} \
         --exclude-filtered

      passing_INDEL="$(gatk CountVariants -V ${basename}.indel.passing.vcf.gz | tail -1)"
      echo "$passing_INDEL" > passing_indel.txt

      #filtered only          
      bcftools isec -C -O v -o ${basename}.indel.filtered.vcf -w1 ${basename}.indel.vcf.gz ${basename}.indel.passing.vcf.gz
      filtered_INDEL="$(gatk CountVariants -V ${basename}.indel.filtered.vcf | tail -1)"
      echo "$filtered_INDEL" > filtered_indel.txt

   }

   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      preemptible: preemptible
   }
   output {

      File output_snp_vcf = "${basename}.snp.vcf.gz"
      File output_snp_vcf_idx = "${basename}.snp.vcf.gz.tbi"

      Int passing_SNP = read_int("passing_snp.txt")
      Int filtered_SNP = read_int("filtered_snp.txt")
      
      Int passing_MNP = read_int("passing_mnp.txt")
      Int filtered_MNP = read_int("filtered_mnp.txt")

      File output_indel_vcf = "${basename}.indel.vcf.gz"
      File output_indel_vcf_idx = "${basename}.indel.vcf.gz.tbi"

      Int passing_INDEL = read_int("passing_indel.txt")
      Int filtered_INDEL = read_int("filtered_indel.txt")
      
   }   
}

#creates an IGV session
#given a list of IGV compatible files (as strings)
#reference is either "hg19" or "hg38"
task GenerateIGVSession {

   input {
      String bloodbiopsydocker
      Array[String] input_files
      String reference_version
      String file_name
      Array[String]? input_names
      Int preemptible
   }

   Array[String] input_names_prefix = if defined(input_names) then prefix('-n ', select_first([input_names])) else []

   command {
      bash /usr/writeIGV.sh ${reference_version} ${sep=" " input_files} ${sep=" " input_names_prefix}  > "${file_name}.xml"
   }
   runtime {
      docker: bloodbiopsydocker
      preemptible: preemptible
   }
   output {
      File igv_session = "${file_name}.xml"
   }
}

task ExtractContam {

   input {
      String bloodbiopsydocker
      File contamination_table
      Int preemptible
   }

   command {
      tail -n1 ~{contamination_table} | cut -f2 > fraction_contamination.txt
   }
   runtime {
      docker: bloodbiopsydocker
      preemptible: preemptible
   }
   output {
      Float fracContam=read_float("fraction_contamination.txt")
   }
}

task VariantsToTable {
   input {
      String bloodbiopsydocker
      File vcf
      File vcf_idx
      String output_name
      File? gatk_override
      Int? memory
      Int disk_pad
      Int preemptible
   }

   Int mem = if defined(memory) then memory * 1000 else 3500
   Int compute_mem = mem - 500
   Int disk_size = ceil(size(vcf, "GB") * 1.5 ) + ceil(size(vcf_idx, "GB")) + disk_pad

   command {

      set -e

      export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

      gatk --java-options "-Xmx${compute_mem}m" VariantsToTable \
         -V ${vcf} \
         -F CHROM -F POS -F ID -F REF -F ALT -F QUAL \
         -F FILTER -ASGF AD -ASGF AF -GF NCount -ASGF SB \
         -SMA \
         -O "${output_name}.table"

   }

   runtime {
      docker: bloodbiopsydocker
      disks: "local-disk " + 500 + " HDD"
      memory: mem + " MB"
      preemptible: preemptible
   }
   output {
      File output_table = "${output_name}.table"
   }
}
