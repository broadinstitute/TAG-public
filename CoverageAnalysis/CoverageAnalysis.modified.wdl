task BaseName {
   File bam_file
   Float bam_size
   Float disk_size
   Int preemptible_tries

   command <<<
      echo ${bam_size} > bam_size.txt

      samtools view -H ${bam_file} | \
      sed -n "/SM:/{s/.*SM:\\(\\)/\\1/; s/\\t.*//p ;q};" | \
      awk -v PREFIX=`basename ${bam_file} ".bam"` '{print PREFIX"_"$0}'
   >>>
   output {
      String base_name_out = read_string(stdout())
      Float sample_bam_size = read_float("bam_size.txt")
   }
   runtime {
      docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      preemptible: preemptible_tries
      memory: "4 GB"
   }
}

task NormalizedCoverage {
  File sample_set
  File coverage_script
  String base_name
  Int preemptible_tries

  Int? extra_disk
  Int disk_size = 50 + select_first([extra_disk,0])

  command {
    Rscript ${coverage_script} ${sample_set} ${base_name}
  }
  output {
    File normalized_coverage = "${base_name}.target_coverage.txt"
    File lc_summary = "${base_name}.lc_targets_summary.txt"
    File lc_list = "${base_name}.lc_targets_list.txt"
    File hist_pdf = "${base_name}.target_coverage_histogram.pdf"
    File hist_png = "${base_name}.target_coverage_histogram.png"
  }
  runtime {
    docker: "us.gcr.io/tag-team-160914/tag-tools:0.3.2"
    disks: "local-disk " + disk_size + " HDD"
    memory: "5GB"
    cpu: "1"
  }
}

task CalculateCoverage {
   File gatk_jar
   File bam_file
   File bam_index
   File interval_list
   File ref_fasta
   File ref_fasta_index
   File ref_dict
   String base_name
   Int min_base_quality
   Int min_mapping_quality
   Float disk_size
   Int preemptible_tries
   String extra_arguments 


   command <<<
      java -jar ${gatk_jar} -T DepthOfCoverage \
           -R ${ref_fasta} \
           -L ${interval_list} \
           -I ${bam_file} \
           -mmq ${min_mapping_quality} \
           -mbq ${min_base_quality} \
           -o "coverage.out" \
           --omitPerSampleStats \
           --printBaseCounts \
           ${extra_arguments}

      # Replace the SM tag with the sample name used in the pipeline
      awk -v BASE_NAME=${base_name} \
          'BEGIN{ FS="\t"; OFS="\t"; }
                { if($1 ~ /^Locus/){$4=BASE_NAME} print $0}' "coverage.out" |
      cut -f1,4 > "${base_name}.coverage"
   >>>
   output {
      File coverage = "${base_name}.coverage"
      Float coverage_file_size = size("${base_name}.coverage", "MB")
   }
   runtime {
      docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      memory: "16 GB"
      preemptible: preemptible_tries
   }
}

task CollectData {
   File picard_jar
   File interval_list
   File sample_genders
   File summary_python_script

   Array[File] coverage_outputs
   Array[File] lcs_pos_outputs
   Array[File] lcs_frac_outputs
   Array[File] mean_coverage_outputs

   String base_name
   String female_chromosome
   String male_chromosome
   Int low_coverage_threshold
   Float sample_threshold
   Float disk_size
   Float? additional_disk_GB

   Float final_disk_size = if defined(additional_disk_GB) then disk_size + additional_disk_GB else disk_size
   Int preemptible_tries

   command <<<
      # Mark samples with unknown genders
      awk '{if($2 == "Unknown"){print($1".");}}' ${sample_genders} > filtered_samples.txt

      FILTER="filtered_samples.txt"

      # Per base coverage matrix
      echo ${sep=' ' coverage_outputs} | tr ' ' '\n' | grep -v -F -f $FILTER > list_
      for F in `cat list_`; do cut -f2 $F > `basename $F`.tmp_; done
      paste <(cut -f1 ${coverage_outputs[0]}) `ls *.tmp_ | sort` > "${base_name}.per_base_coverage.txt" && rm *.tmp_

      # Target mean coverage
      echo ${sep=' ' mean_coverage_outputs} | tr ' ' '\n' | grep -v -F -f $FILTER > list_
      for F in `cat list_`; do cut -f2 $F > `basename $F`.tmp_; done
      paste <(cut -f1 ${mean_coverage_outputs[0]}) `ls *.tmp_ | sort` > "${base_name}.mean_coverage.txt" && rm *.tmp_

      # Target low covered base fraction
      echo ${sep=' ' lcs_frac_outputs} | tr ' ' '\n' | grep -v -F -f $FILTER > list_
      for F in `cat list_`; do cut -f2 $F > `basename $F`.tmp_; done
      paste <(cut -f1 ${lcs_frac_outputs[0]}) `ls *.tmp_ | sort` > "${base_name}.low_covered_frac.txt" && rm *.tmp_

      # Per base low covered site matrix
      echo ${sep=' ' lcs_pos_outputs} | tr ' ' '\n' | grep -v -F -f $FILTER > list_
      for F in `cat list_`; do cut -f3 $F > `basename $F`.tmp_; done
      paste <(cut -f1,2 ${lcs_pos_outputs[0]}) `ls *.tmp_ | sort` > low_covered_flags.txt && rm *.tmp_

      # Filtered samples with unknown gender
      if [[ `cat $FILTER | wc -l` -eq 0 ]]
      then
         echo "None." > ${base_name}.filtered_samples.txt
         touch ${base_name}.filtered_samples.per_base_coverage.txt
      else
         echo ${sep=' ' coverage_outputs} | tr ' ' '\n' | grep -F -f $FILTER > list_
         for F in `cat list_`; do cut -f2 $F > `basename $F`.tmp_; done
         paste <(cut -f1 ${coverage_outputs[0]}) `ls *.tmp_ | sort` > "${base_name}.filtered_samples.per_base_coverage.txt" && rm *.tmp_
         sed 's/\.$//g' $FILTER > ${base_name}.filtered_samples.txt
      fi

      python ${summary_python_script} \
             ${sample_threshold} \
             ${base_name}.per_base_coverage.txt \
             low_covered_flags.txt > ${base_name}.sample_set.txt
   >>>
   output {
      File sample_set_summary = "${base_name}.sample_set.txt"
      File coverage_matrix = "${base_name}.per_base_coverage.txt"
      File mean_coverage_matrix = "${base_name}.mean_coverage.txt"
      File lcs_fraction_matrix = "${base_name}.low_covered_frac.txt"
      File filtered_samples = "${base_name}.filtered_samples.txt"
      File filtered_samples_coverage = "${base_name}.filtered_samples.per_base_coverage.txt"
   }
   runtime {
      docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
      disks: "local-disk " + sub(final_disk_size, "\\..*", "") + " HDD"
      memory: "64 GB"
      preemptible: preemptible_tries
   }
}

task SumFloats {
   Array[Float] sizes
   Int preemptible_tries

   command <<<
      python -c "print ${sep='+' sizes}"
   >>>
   output {
      Float total_size = read_float(stdout())
   }
   runtime {
     disks: "local-disk " + 10 + " HDD"
     memory: "5000 MB"
     docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
     preemptible: preemptible_tries
   }
}

task CollectGenderWES {
   File picard_jar
   Array[File] bam_files
   Array[File] bam_indices
   Array[String] base_names
   String male_chromosome
   String female_chromosome
   String output_prefix
   Float disk_size
   Int? extra_disk
   Float total_disk = disk_size + select_first([extra_disk,0])
   Int preemptible_tries

   command <<<

      # Obtain file names for Picard input
      echo ${sep=' ' bam_files} | tr ' ' '\n' > bam.list

      # Extract SM tags and combine these with sample base names
      paste <(echo ${sep=' ' base_names} | tr ' ' '\n') <(cat bam.list | xargs -I % sh -c 'basename % .bam') | \
      xargs -n2 python -c 'import sys; print sys.argv[1].split(sys.argv[2]+"_")[1]' | \
      paste - <(echo ${sep=' ' base_names} | tr ' ' '\n') > sample_names.txt

      # Run Picard tool to estimate genders
      java -jar ${picard_jar} InferSexFromBAM \
                              MALE_CHROMS=${male_chromosome} \
                              FEMALE_CHROMS=${female_chromosome} \
                              I=bam.list O=/dev/stdout | cut -f1,5 > picard.genders.out

      # Python snippet to create a gender table with sample base names
      ## ==== python code block starts ====
PYTHON_SRC="
import sys
_gender_out = sys.argv[1]
_sample_tab = sys.argv[2]
gender = {'1':'Male', '2':'Female'}
estimation = dict([ x.rstrip().split('\t') for x in open(_gender_out) ])
for x in open(_sample_tab):
    x = x.rstrip().split('\t')
    print(x[1] + '\t' + gender[estimation[x[0]]])
"
      ## ==== python code block ends ====

      # Create a final gender table
      python -c "$PYTHON_SRC" picard.genders.out sample_names.txt > ${output_prefix}.genders.txt
   >>>
   output {
      File gender_out = "${output_prefix}.genders.txt"
   }
   runtime {
      docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
      disks: "local-disk " + sub(total_disk, "\\..*", "") + " HDD"
      memory: "30 GB"
      preemptible: preemptible_tries
   }
}

task DetermineGender {
   File determine_gender_r_script
   File bam_file
   File bam_index
   File autosome_list
   String male_chromosome
   String female_chromosome
   String base_name
   Float disk_size
   Int preemptible_tries

   command <<<
      samtools idxstats ${bam_file} > ${base_name}.idxstats.txt
      Rscript -e "source('${determine_gender_r_script}'); \
                  determineGender('${base_name}.idxstats.txt', '${base_name}', '${base_name}.genders.txt', \
                                  '${autosome_list}', '${male_chromosome}', '${female_chromosome}')"
   >>>
   output {
      File gender_out = "${base_name}.genders.txt"
   }
   runtime {
      docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      memory: "16 GB"
      preemptible: preemptible_tries
   }
}

task CollectGenderWGS {
   Array[File?] genders
   String output_prefix
   Int preemptible_tries

   command <<<
      set -e
      set -v

      # Put all the genders into a single file
      cat ${sep=" " genders} > ${output_prefix}.genders.txt
   >>>
   output {
      File gender_out = "${output_prefix}.genders.txt"
   }
   runtime {
      docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
      disks: "local-disk 10 HDD"
      memory: "4 GB"
      preemptible: preemptible_tries
   }
}

task SampleCoverageStats {
   File picard_jar
   String base_name
   File sample_coverage
   File gender_estimation
   File interval_list
   String male_chromosome
   String female_chromosome
   Int low_coverage_threshold
   Float disk_size
   Int preemptible_tries

   command <<<
      # Convert interval_list to BED
      java -jar ${picard_jar} IntervalListToBed I=${interval_list} O=/dev/stdout | \
      cut -f1-4 > intervals.bed

      # Intersect coverage output with the interval list
      # Output: [1]chromosome [2]start [3]target_name [4]coverage
      grep -v "^Locus" ${sample_coverage} | \
      awk 'BEGIN{ OFS="\t"; FS="\t"; }{ split($1,pos,":"); $1=pos[1]"\t"pos[2]-1"\t"pos[2]; print($0); }' | \
      bedtools intersect -a intervals.bed -b stdin -wo | \
      awk 'BEGIN{ OFS="\t"; FS="\t"; }{ print($5,$7,$4,$8); }' > intersected.txt

      # Extract sample gender
      GENDER=`grep -P '${base_name}\t' ${gender_estimation} | cut -f2`

      # Flag low covered sites
      # Output: [1]chromosome:start [2]target_name [3]low_covered (yes=1 / no=0)
      awk -v gender=$GENDER \
          -v male_chrom="${male_chromosome}" \
          -v female_chrom="${female_chromosome}" \
          -v cutoff=${low_coverage_threshold} \
          'BEGIN{ OFS="\t"; FS="\t"; }
          {  c = cutoff;
              if(gender == "Male" && ($1 == female_chrom || $1 == male_chrom)){
                 c = cutoff/2; # Use a half of the cutoff
              }
              res = (c > $4); # Mark low covered sites
              if(gender == "Female" && $1 == male_chrom){
                 res = 1; # Mark sites on male chromosome in female sample as low covered
              }
              print($1":"$2,$3,res);
           }' intersected.txt | \
      sort -k1,1 -k2,2n > ${base_name}.lcs_pos.txt

      # Compute mean coverage
      awk '{ target_len[$3] += 1; sum_coverage[$3] += $4; }
           END{ for(i in target_len){
                   printf("%s\t%.4f\n", i, sum_coverage[i]/target_len[i]);
                }
           }' intersected.txt | \
      sort -k1,1 | cat <(echo -e "Target\t${base_name}") - > ${base_name}.mean_cov.txt

      # Compute low covered fraction
      awk '{ target_len[$2] += 1; lcs[$2] += $3; }
          END{ for(i in target_len){
                   printf("%s\t%.8f\n", i, lcs[i]/target_len[i]);
               }
          }' ${base_name}.lcs_pos.txt | \
      sort -k1,1 | cat <(echo -e "Target\t${base_name}") - > ${base_name}.lcs_frac.txt
   >>>
   output {
      File sample_lcs_position = "${base_name}.lcs_pos.txt"
      File sample_lcs_fraction = "${base_name}.lcs_frac.txt"
      File sample_mean_coverage = "${base_name}.mean_cov.txt"
   }
   runtime {
      docker: "us.gcr.io/tag-team-160914/tag-tools:0.0.1"
      disks: "local-disk " + sub(disk_size, "\\..*", "") + " HDD"
      memory: "16 GB"
      preemptible: preemptible_tries
   }
}

task PlotData {
   File plotting_r_script
   File sample_genders
   File mean_coverage_matrix
   File lcs_fraction_matrix
   File sample_set_summary

   Float? memory_size

   String base_name
   String male_chromosome
   String female_chromosome
   String mitochondrial_dna

   Int low_coverage_threshold
   Int preemptible_tries

   command <<<
      Rscript --vanilla ${plotting_r_script} ${low_coverage_threshold} ${base_name} \
                                            '${male_chromosome}' '${female_chromosome}' '${mitochondrial_dna}' \
                                             ${sample_genders} ${mean_coverage_matrix} ${lcs_fraction_matrix} \
                                             ${sample_set_summary} || touch ${base_name}.callable_frac.pdf ${base_name}.target_coverage.pdf
   >>>
   output {
      File callable_frac_pdf = "${base_name}.callable_frac.pdf"
      File target_coverage_pdf = "${base_name}.target_coverage.pdf"
   }
   runtime {
      docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
      disks: "local-disk 10 HDD"
      memory: select_first([memory_size, 4]) + " GB"
      preemptible: preemptible_tries
   }
}

task calcMeanCov {
	Array[File] mean_coverage_outputs
    File sample_genders
    
    command <<<
    	awk '{if($2 == "Unknown"){print($1".");}}' ${sample_genders} > filtered_samples.txt

      	FILTER="filtered_samples.txt"
    
    	echo ${sep=' ' mean_coverage_outputs} | tr ' ' '\n' | grep -v -F -f $FILTER > list_
      	for F in `cat list_`; do cut -f2 $F > `basename $F`.tmp_; done
      	paste <(cut -f1 ${mean_coverage_outputs[0]}) `ls *.tmp_ | sort` > "sample_mean_coverage.txt" && rm *.tmp_
      
    	cut -f2- sample_mean_coverage.txt | tail -n+2 | \
        awk -v sum=0 -v count=0 'BEGIN{FS="\t"}{for(i=1; i<=NF; i++){sum+=$i; count+=1}}END{print(sum/count)}' > mean_coverage.txt
    >>>
    
    output{
    	Float mean_cov = read_float("mean_coverage.txt")
    }
    
    runtime{
      docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
      disks: "local-disk 10 HDD"
      memory: "2.5 GB"
      preemptible: 3
    }
}

workflow ScatteredCoverage {
   File gatk_jar
   File picard_jar

   Boolean is_exome
   Boolean plot_coverage
   File summary_python_script
   File determine_gender_r_script
   File plotting_r_script
   File autosome_list
   String male_chromosome
   String female_chromosome
   String mitochondrial_dna

   String output_prefix
   Array[File] bam_files
   Array[File] bam_indices
   File interval_list
   File ref_fasta
   File ref_fasta_index
   File ref_dict
   Int min_base_quality
   Int min_mapping_quality
   Int low_coverage_threshold
   Float sample_fraction_threshold

   Int preemptible_tries
   Int? increase_disk_size
   Int additional_disk = select_first([increase_disk_size, 20])

   Array[Pair[File, File]] bam_index_pairs = zip(bam_files, bam_indices)
   Float sample_threshold = sample_fraction_threshold * length(bam_files)

   File coverage_script
   Boolean run_normalized_coverage

   scatter (paired_data in bam_index_pairs) {

      Float bam_size = size(paired_data.left, "GB")

      call BaseName {
         input: bam_file = paired_data.left,
                bam_size = bam_size,
                disk_size = bam_size + additional_disk,
                preemptible_tries = preemptible_tries
      }

      call CalculateCoverage {
         input: gatk_jar = gatk_jar,
                bam_file = paired_data.left,
                bam_index = paired_data.right,
                base_name = BaseName.base_name_out,
                interval_list = interval_list,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                min_base_quality = min_base_quality,
                min_mapping_quality = min_mapping_quality,
                disk_size = (bam_size * 1.5) + additional_disk,
                preemptible_tries = preemptible_tries,
                extra_arguments = "--countType COUNT_FRAGMENTS_REQUIRE_SAME_BASE -allowPotentiallyMisencodedQuals"
                
      }

      if (!is_exome) {
         call DetermineGender {
            input: determine_gender_r_script = determine_gender_r_script,
                   bam_file = paired_data.left,
                   bam_index = paired_data.right,
                   autosome_list = autosome_list,
                   male_chromosome = male_chromosome,
                   female_chromosome = female_chromosome,
                   base_name = BaseName.base_name_out,
                   disk_size = bam_size + additional_disk,
                   preemptible_tries = preemptible_tries
         }
      }
   }

   if (!is_exome) {
      call CollectGenderWGS {
         input: genders = DetermineGender.gender_out,
                output_prefix = output_prefix,
                preemptible_tries = preemptible_tries
      }
   }

   if (is_exome) {
      call SumFloats as ExomeBamSize {
         input: sizes = BaseName.sample_bam_size,
                preemptible_tries = preemptible_tries
      }
      call CollectGenderWES {
         input: picard_jar = picard_jar,
                bam_files = bam_files,
                bam_indices = bam_indices,
                base_names = BaseName.base_name_out,
                male_chromosome = male_chromosome,
                female_chromosome = female_chromosome,
                output_prefix = output_prefix,
                disk_size = ExomeBamSize.total_size + additional_disk,
                preemptible_tries = preemptible_tries
      }
   }

   Array[Pair[File,Float]] cvr = zip(CalculateCoverage.coverage, CalculateCoverage.coverage_file_size)
   Array[Pair[String,Pair[File,Float]]] sample_paired_data = zip(BaseName.base_name_out, cvr)

   scatter (sample_data in sample_paired_data) {

      Pair[File, Float] coverage_data = sample_data.right

      call SampleCoverageStats {
         input: picard_jar = picard_jar,
                base_name = sample_data.left,
                sample_coverage = coverage_data.left,
                gender_estimation = gender_estimation,
                interval_list = interval_list,
                male_chromosome = male_chromosome,
                female_chromosome = female_chromosome,
                low_coverage_threshold = low_coverage_threshold,
                disk_size = (coverage_data.right/1000) + additional_disk,
                preemptible_tries = preemptible_tries
      }
   }
   if(run_normalized_coverage){
   	 call calcMeanCov{
       input: mean_coverage_outputs = SampleCoverageStats.sample_mean_coverage,
       		  sample_genders = gender_estimation
     }
     
     Int baitQC_threshold = ceil(0.2*calcMeanCov.mean_cov)
     
     scatter (sample_data in sample_paired_data) {

      	Pair[File, Float] coverage_data_updateLC = sample_data.right

      	call SampleCoverageStats as SampleCoverageStats_updateLC {
         	input: picard_jar = picard_jar,
                   base_name = sample_data.left,
                   sample_coverage = coverage_data_updateLC.left,
                   gender_estimation = gender_estimation,
                   interval_list = interval_list,
                   male_chromosome = male_chromosome,
                   female_chromosome = female_chromosome,
                   low_coverage_threshold = baitQC_threshold,
                   disk_size = (coverage_data_updateLC.right/1000) + additional_disk,
                   preemptible_tries = preemptible_tries
        }
   	 }
     
   }

   call SumFloats as CoverageFileSize {
      input: sizes = CalculateCoverage.coverage_file_size,
             preemptible_tries = preemptible_tries
   }

   File gender_estimation = select_first([CollectGenderWES.gender_out, CollectGenderWGS.gender_out])

   call CollectData {
      input: picard_jar = picard_jar,
             interval_list = interval_list,
             sample_genders = gender_estimation,
             summary_python_script = summary_python_script,
             coverage_outputs = CalculateCoverage.coverage,
             lcs_pos_outputs = select_first([SampleCoverageStats_updateLC.sample_lcs_position,SampleCoverageStats.sample_lcs_position]),
             lcs_frac_outputs = select_first([SampleCoverageStats_updateLC.sample_lcs_fraction,SampleCoverageStats.sample_lcs_fraction]),
             mean_coverage_outputs = SampleCoverageStats.sample_mean_coverage,
             base_name = output_prefix,
             female_chromosome = female_chromosome,
             male_chromosome = male_chromosome,
             low_coverage_threshold = select_first([baitQC_threshold,low_coverage_threshold]),
             sample_threshold = sample_threshold,
             disk_size = (CoverageFileSize.total_size/1000) + additional_disk,
             preemptible_tries = preemptible_tries
   }
   if(run_normalized_coverage){
   	 
     call NormalizedCoverage{
       input: sample_set = CollectData.sample_set_summary,
              base_name = output_prefix,
              preemptible_tries = preemptible_tries,
              coverage_script = coverage_script
     }
   }

  if(plot_coverage){
  		call PlotData{
      		input: plotting_r_script = plotting_r_script,
            	   sample_genders = gender_estimation,
             	   mean_coverage_matrix = CollectData.mean_coverage_matrix,
             	   lcs_fraction_matrix = CollectData.lcs_fraction_matrix,
             	   sample_set_summary = CollectData.sample_set_summary,
             	   base_name = output_prefix,
             	   male_chromosome = male_chromosome,
             	   female_chromosome = female_chromosome,
             	   mitochondrial_dna = mitochondrial_dna,
             	   low_coverage_threshold = select_first([baitQC_threshold,low_coverage_threshold]),
             	   preemptible_tries = preemptible_tries
   		}
	}


   output {
      File coverage_output = CollectData.coverage_matrix
      File sample_genders = "${gender_estimation}"
      File sample_set_summary = CollectData.sample_set_summary
      File mean_coverage_matrix = CollectData.mean_coverage_matrix
      File lcs_fraction_matrix = CollectData.lcs_fraction_matrix
      File filtered_samples = CollectData.filtered_samples
      File filtered_sample_coverage = CollectData.filtered_samples_coverage
      
      File? callable_frac_pdf = PlotData.callable_frac_pdf
      File? target_coverage_pdf = PlotData.target_coverage_pdf
      
      File? normalized_coverage_table = NormalizedCoverage.normalized_coverage
      File? coverage_histogram = NormalizedCoverage.hist_png
      File? lc_sample_summary = NormalizedCoverage.lc_summary
      File? lc_sample_list = NormalizedCoverage.lc_list
   }
}
