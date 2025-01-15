version 1.0

workflow Callability_Analysis {
    input {
        Array[File] bam_or_cram_paths
        Array[File] bai_or_crai_paths
        Array[File] sample_ids
        Int min_base_quality
        Int min_mapping_quality
        Int low_coverage_threshold
        Int sample_fraction
        File target_bed
        File ref_dict
        File ref_fasta
        File ref_index
        File ref_pac
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_sa
        Int coverage_mem_gb
        String chrX_name
        String chrY_name
        String output_prefix
        File gene_bed
        File exome_bed
        File? gene_list
        Boolean generate_gene_summary
        Boolean generate_plot = false
    }
    
   Float one = 1.0

   Array[Pair[File, File]] bam_index_pairs = zip(bam_or_cram_paths, bai_or_crai_paths)

   scatter (name_and_paired_data in zip(sample_ids, bam_index_pairs)) {
      String sample_name = name_and_paired_data.left
	  Pair[File, File] paired_data = name_and_paired_data.right

      call DetermineXYCoverage {
         input:
            bam_or_cram_path = paired_data.left,
            bai_or_crai_path = paired_data.right,
            chrX_name = chrX_name,
            chrY_name = chrY_name,
            base_name = sample_name,
            disk_size = 200
      }
      call CalculateCoverage {
          input:
            bam_or_cram_path = paired_data.left,
            bai_or_crai_path = paired_data.right,
            sample_name = sample_name,
            target_bed = target_bed,
            min_base_quality = min_base_quality,
            min_mapping_quality = min_mapping_quality,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_index = ref_index,
            coverage_mem_gb = coverage_mem_gb
      }
   }


   call DetermineSexesByClustering {
        input:
            inputFiles = DetermineXYCoverage.xy_out
   }

   call CollectData {
      input:
        coverage_files = CalculateCoverage.coverage,
        sexes_file = DetermineSexesByClustering.output_sexes,
        chrX_name = chrX_name,
        chrY_name = chrY_name,
        base_name = output_prefix,
        low_coverage_threshold = low_coverage_threshold,
        sample_threshold = sample_fraction/100.0
   }
   call BedToIntervalList as InputBedToInterval {
        input:
            bed = target_bed,
            refDict = ref_dict

   }

   call BedToIntervalList as UndercoveredBedToInterval {
        input:
            bed = CollectData.bed_file_output,
            refDict = ref_dict
    }
    
   call CountBases as CountUndercoveredBases {
        input:
            intervalListOrVcf = UndercoveredBedToInterval.interval_list
    }

   call CountBases as CountAllBases {
        input:
            intervalListOrVcf = InputBedToInterval.interval_list
    }
        
    File undercovered_bed = CollectData.bed_file_output
    File coverage_output = CollectData.coverage_output
    File undercovered_interval_list = UndercoveredBedToInterval.interval_list
    Int bases_undercovered = CountUndercoveredBases.bases
    Float fraction_undercovered =  CountUndercoveredBases.bases/(CountAllBases.bases/one)
    
    call CheckFileNotEmpty {
        input:
            file_to_check = CollectData.bed_file_output
    }
    
    if (CheckFileNotEmpty.is_not_empty) {
        call GenerateAnnotation {
            input:
            exome_bed = exome_bed,
            bed_to_annotate = CollectData.bed_file_output,
            output_prefix = output_prefix,
            gene_list = gene_list,
            gene_bed = gene_bed
        }

        scatter (i in range(length(sample_ids))) {
            call samtools_coverage {
                input:
                    bam_or_cram_path = bam_or_cram_paths[i],
                    bai_or_crai_path = bai_or_crai_paths[i],
                    bed_to_annotate = CollectData.bed_file_output,
                    sample_id = sample_ids[i],
                    ref_dict = ref_dict,
                    ref_fasta = ref_fasta,
                    ref_index = ref_index,
                    ref_pac = ref_pac,
                    ref_amb = ref_amb,
                    ref_ann = ref_ann,
                    ref_bwt = ref_bwt,
                    ref_sa = ref_sa
            }
        }

        call summarize_coverage {
            input:
            coverage_outputs = samtools_coverage.coverage_output,
            sample_fraction = sample_fraction
        }

        call generate_clinvar_results {
            input:
            bed_to_annotate = CollectData.bed_file_output
        }

        call concatenate_results {
            input:
            clinvar_annotation = generate_clinvar_results.clinvar_annotation, 
            samtools_coverage_file = summarize_coverage.samtools_coverage_summary,
            annotation_file = GenerateAnnotation.annotation_per_interval,
            output_prefix = output_prefix
        }
    if (generate_gene_summary) {
        call GenerateGeneSummary {
            input:
            CoverageFile = CollectData.coverage_output,
            gene_bed = gene_bed,
            group_by_gene = GenerateAnnotation.grouped_by_gene,
            sample_fraction = sample_fraction,
            generate_plot = generate_plot
     }
     }
  }

output {
        File out_undercovered_bed = undercovered_bed
        File out_undercovered_interval_list = undercovered_interval_list
        File out_coverage_output = coverage_output
        Int out_bases_undercovered = bases_undercovered
        Float out_fraction_undercovered = fraction_undercovered
        File? annotation_per_interval = GenerateAnnotation.annotation_per_interval
        File? grouped_by_gene = GenerateAnnotation.grouped_by_gene
        Array[File]? coverage_outputs = samtools_coverage.coverage_output
        File? samtools_coverage_summary = summarize_coverage.samtools_coverage_summary
        File? clinvar_annotation = generate_clinvar_results.clinvar_annotation
        File? integrated_annotation_file = concatenate_results.integrated_annotation_file
        File? callable_fraction_plot = GenerateGeneSummary.output_plot
        File? uncallable_gene_output_table = GenerateGeneSummary.gene_base_count
    }
}

task CalculateCoverage {
    input {
       File bam_or_cram_path
       File bai_or_crai_path
       File target_bed
       Int min_base_quality
       Int min_mapping_quality
       Int disk_size = 200
       String sample_name
       File ref_fasta
       File ref_dict
       File ref_index
       Int coverage_mem_gb
   }

   command <<<
      set -e

    echo ~{bam_or_cram_path} > alignment_file.txt
      
    samtools depth -q ~{min_base_quality} -Q ~{min_mapping_quality} -b ~{target_bed} -aa -H -J -s -f alignment_file.txt -o ~{sample_name}.coverage
    awk 'BEGIN{FS="\t"; OFS="\t"}
        NR==1 {
            print "Locus", "~{sample_name}"; 
            next
        }
        {
            printf $1":"$2;  # Combine first two columns to form Locus
            split($3, a, ".");
            printf "\t"a[1]; 
            printf "\n";
        }' ~{sample_name}.coverage > ~{sample_name}.coverage.txt
      
   >>>
   output {
      File coverage = "~{sample_name}.coverage.txt"
   }
   runtime {
      docker: "us.gcr.io/tag-public/samtools_r:v1"
      disks: "local-disk " + disk_size + " HDD"
      memory: "~{coverage_mem_gb} GB"
   }
}

task CollectData {
	input {
       Array[File] coverage_files
       String chrX_name
       String chrY_name
       File sexes_file
       String base_name
       Int low_coverage_threshold
       Float sample_threshold
       Int disk_size = 200
       Int memory_gb = 32
   }

   command <<<
      set -e
      set -v
      
      echo "Merging all coverage files into one"
      
      COVERAGE_FILES="~{sep=' ' coverage_files}"

      first_file=$(echo ${COVERAGE_FILES} | awk '{print $1}')

      cp "$first_file" ~{base_name}.coverage.txt

      for file in $(echo ${COVERAGE_FILES} | cut -d' ' -f2-); do

          paste ~{base_name}.coverage.txt <(awk 'BEGIN {FS=OFS="\t"} {print $2}' "$file") > temp.txt
          mv temp.txt ~{base_name}.coverage.txt
      done
        
      echo "Starting Data Munging in R"

      # In R do various data munging operations
      Rscript -<<"EOF"
      library(data.table)
      options(scipen = 999)
      sites <- fread("~{base_name}.coverage.txt")
      sexes <- fread("~{sexes_file}", col.names = c("sample", "sex", "normalizedX", "normalizedY"))
      maleColumns <- sexes[sex == "Male", sample]
      femaleColumns <- sexes[sex == "Female", sample]

      sites[, contig := gsub(":.*", "", Locus)]
      sites[, pos := as.numeric(gsub(".*:", "", Locus))]
      sites[, countBelowThreshold := Reduce(`+`, lapply(.SD, "<", ~{low_coverage_threshold})), .SDcols = c(maleColumns, femaleColumns)]
      sites[contig == "~{chrX_name}", maleCountBelowThreshold := Reduce(`+`, lapply(.SD, "<", ~{low_coverage_threshold} / 2)), .SDcols = maleColumns]
      sites[contig == "~{chrX_name}", femaleCountBelowThreshold := Reduce(`+`, lapply(.SD, "<", ~{low_coverage_threshold})), .SDcols = femaleColumns]
      sites[contig == "~{chrX_name}", countBelowThreshold := femaleCountBelowThreshold + maleCountBelowThreshold]
      sites[contig == "~{chrY_name}", countBelowThreshold := Reduce(`+`, lapply(.SD, "<", ~{low_coverage_threshold} / 2)), .SDcols = maleColumns]

      sampleThreshold <- sexes[, .N] * ~{sample_threshold}
      maleSampleThreshold <- sexes[sex == "Male", .N] * ~{sample_threshold}

      sites[contig != "~{chrY_name}", poorlyCovered := countBelowThreshold >= sampleThreshold]
      sites[contig == "~{chrY_name}", poorlyCovered := countBelowThreshold >= maleSampleThreshold]
      poorlyCoveredSites <- sites[poorlyCovered == TRUE, .(contig, pos)]
      poorlyCoveredSites[, interval := cumsum(c(1, (pos - shift(pos) > 1 | contig != shift(contig))[-1]))]
      poorlyCoveredIntervals <- unique(poorlyCoveredSites[, .(contig = contig, start = min(pos) - 1, end = max(pos)), by = interval])
      poorlyCoveredIntervals[, orientation := "+"]
      poorlyCoveredIntervals[, name := paste0(contig, "_", start, "_", end)]
      poorlyCoveredIntervals[, interval := NULL]
      poorlyCoveredIntervals
      fwrite(poorlyCoveredIntervals, "~{base_name}.low_coverage.bed", col.names = FALSE, sep = "\t")
      EOF

   >>>
   output {
      File coverage_output = "~{base_name}.coverage.txt"
      File bed_file_output = "~{base_name}.low_coverage.bed"
   }
   runtime {
      docker: "us.gcr.io/tag-public/samtools_r:v1"
      disks: "local-disk " + disk_size + " HDD"
      memory: memory_gb + "GB"
   }
}

task DetermineXYCoverage {
	input {
        File bam_or_cram_path
        File bai_or_crai_path
        String chrX_name
        String chrY_name
        String base_name
        Int disk_size = 200
      }

    command <<<
         samtools idxstats ~{bam_or_cram_path} > ~{base_name}.idxstats.txt

         Rscript -<<"EOF"
         idxstats <- read.table("~{base_name}.idxstats.txt", col.names = c("contig", "length", "mapped", "unmapped"))

          # Determine normalization constant.  This is the number of bases/chromosome/read.
         normalization <- sum(as.numeric(idxstats[grepl("^chr[0-9]{1}$|^chr[0-9]{2}$", idxstats$contig), ]$length)) /
              sum(as.numeric(idxstats[grepl("^chr[0-9]{1}$|^chr[0-9]{2}$", idxstats$contig), ]$mapped))
          # Estimate the number of X and Y chromosomes in the sample by
          # normalizing the coverage on X and Y.
         normalizedX <- 2 / (idxstats[idxstats$contig == "~{chrX_name}", ]$length /
                                  idxstats[idxstats$contig == "~{chrX_name}", ]$mapped / normalization)
         normalizedY <- 2 / (idxstats[idxstats$contig == "~{chrY_name}", ]$length /
                                  idxstats[idxstats$contig == "~{chrY_name}", ]$mapped / normalization)
         sampleAndXYChroms <- data.frame("~{base_name}", round(normalizedX, 3), round(normalizedY, 3))

         write.table(sampleAndXYChroms, "~{base_name}.xy.txt", col.name = F, row.name = F, quote = F)
         EOF
    >>>

    output {
          File xy_out = "~{base_name}.xy.txt"
          File idxstats = "~{base_name}.idxstats.txt"
    }

    runtime {
          docker: "us.gcr.io/tag-public/samtools_r:v1"
          disks: "local-disk " + disk_size + " HDD"
          memory: "16 GB"
       }
}

task DetermineSexesByClustering {
	input {
    	Array[File] inputFiles
        Int disk_size = 32
    }

    command <<<
    Rscript -<<"EOF"
    library(data.table)

    d_list <- lapply(list("~{sep="\",\"" inputFiles}"), function(file) fread(file, colClasses = c("V1" = "character")))
    d <- rbindlist(d_list)
    d[, V1 := as.character(V1)]
    d_mat <- as.matrix(d,rownames = "V1")
    initial_centers <- matrix(c(1,2,1,0), nrow=2)
    clusters <- kmeans(d_mat,initial_centers)
    d[,cluster:=list(clusters$cluster)]
    d_centers <- data.table(clusters$centers)
    d_centers[,d_male:=sqrt((V2-1)^2 + (V3-1)^2)]
    d_centers[,d_female:=sqrt((V2-2)^2 + (V3-0)^2)]
    if (d_centers[1,d_male]<d_centers[2,d_male] && d_centers[1,d_female]>d_centers[2,d_female]) {
        d[cluster==2,sex:="Female"]
        d[cluster==1,sex:="Male"]
    } else if (d_centers[1,d_male]>d_centers[2,d_male] && d_centers[1,d_female]<d_centers[2,d_female]) {
        d[cluster==2,sex:="Male"]
        d[cluster==1,sex:="Female"]
    } else {
        d[,sex:="Unknown"]
    }
    d[,cluster:=NULL]
    setcolorder(d,c("V1","sex","V2","V3"))
    write.table(d, "all.sexes.txt", col.name = F, row.name = F, quote = F)
    EOF

    >>>

    runtime {
        docker: "us.gcr.io/tag-public/samtools_r:v1"
        disks: "local-disk " + disk_size + " HDD"
        memory: "16 GB"
    }

    output {
        File output_sexes = "all.sexes.txt"
    }
}

task BedToIntervalList {
	input {
        File bed
        File refDict
    }

    String intervalListOut = sub(basename(bed), ".bed$", ".interval_list")

    Int disk_size = 10 + 2*ceil(size(bed, "GB") + size(refDict, "GB"))
    command <<<
        java -jar /usr/gitc/picard.jar BedToIntervalList I=~{bed} O=~{intervalListOut} SD=~{refDict}
    >>>

    runtime {
      docker: "broadinstitute/genomes-in-the-cloud:2.2.5-1486412288"
      memory: "4 GB"
      disks: "local-disk " + disk_size + " HDD"
   }

   output {
   File interval_list = "${intervalListOut}"
   }
}

task CountBases {
    input {
        File intervalListOrVcf
    }

    Int disk_size = 10 + ceil(size(intervalListOrVcf, "GB"))
    File picardJar = "gs://gptag/AnnotateBed/picard.jar"

    command <<<
        if [[ ~{intervalListOrVcf} == *vcf ]]; then
            java -jar /usr/gitc/picard.jar VcfToIntervalList I=~{intervalListOrVcf} O=vcf.interval_list
            java -jar ~{picardJar} IntervalListTools I=vcf.interval_list COUNT_OUTPUT=bases.txt OUTPUT_VALUE=BASES
        else
            java -jar ~{picardJar} IntervalListTools I=~{intervalListOrVcf} COUNT_OUTPUT=bases.txt OUTPUT_VALUE=BASES
        fi
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
        preemptible: 0
        disks: "local-disk " + disk_size + " HDD"
        memory: "16 GB"
    }

    output {
        Int bases = read_int("bases.txt")
    }
}

task CheckFileNotEmpty {
    input {
      File file_to_check
    }

    command {
        # Check if file is not empty
        if [ -s "~{file_to_check}" ]; then
            echo true
        else
            echo false
        fi
    }
    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
        preemptible: 0
        disks:  "local-disk 4 HDD"
        memory: "4 GB"
    }

    output {
        Boolean is_not_empty = read_boolean(stdout())
    }
}

task GenerateAnnotation {
    input {
        File bed_to_annotate
        File gene_bed
        File exome_bed
        File? gene_list
        String output_prefix
        Int maxRetries = 1
        Int preemptible = 3
        Int diskGB = 50
        String? docker_override
    }
    command {
        python3 /scripts/bed_annotate_script_hg38.py --annotation /reference_files/MANE.GRCh38.v1.3.ensembl_genomic.gtf --bed ~{bed_to_annotate} --gene_bed ~{gene_bed} --exome_bed ~{exome_bed} ~{'--gene_list ' + gene_list} --output_prefix ~{output_prefix}
    }
    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/annotatebed_hg38:v3"])
        preemptible: preemptible
        disks: "local-disk ~{diskGB} HDD"
        memory: "8 GB"
        maxRetries: maxRetries
    }
    output {
        File annotation_per_interval = "~{output_prefix}.grouped_by_interval.annotated.txt"
        File grouped_by_gene = "~{output_prefix}.grouped_by_gene.txt"
    }
}

task samtools_coverage {
    input {
      File bam_or_cram_path
      File bai_or_crai_path
      File bed_to_annotate
      File ref_fasta
      File ref_index
      File ref_dict
      File ref_pac
      File ref_amb
      File ref_ann
      File ref_bwt
      File ref_sa
      String sample_id
      Int memory_gb = 64
      Int disk_size = 100
    }
    command <<<
        for region in `awk '{print $1":"$2"-"$3}' ~{bed_to_annotate}`
        do
           samtools coverage ~{bam_or_cram_path} -r ${region} --reference ~{ref_fasta} >> tmp
        done
        awk '/startpos/&&c++>0 {next} 1' tmp > ~{sample_id}.coverage.txt
>>>

    output {
        File coverage_output = "~{sample_id}.coverage.txt"
    }
    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
        memory: memory_gb + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
task summarize_coverage {
    input{
        Array[File] coverage_outputs
        Int memory_gb = 16
        Int disk_size = 32
        Int sample_fraction
        String? docker_override
    }
    command {
    python3 /scripts/summarize_coverage.hg38.py "~{sep=' ' coverage_outputs}" ~{sample_fraction}
    }
    output{
    File samtools_coverage_summary = 'samtools_coverage_summary.txt'
    }
    runtime {
    	docker: select_first([docker_override, "us.gcr.io/tag-public/annotatebed_hg38:v3"])
        memory: memory_gb + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
task generate_clinvar_results{
    input {
        File bed_to_annotate
        Int memory_gb = 16
        Int disk_size = 32
        String? docker_override
    }
    command <<<
        sed 's/^chr//' ~{bed_to_annotate} > non_chr_bed.tmp
        bcftools query -f '%CHROM\t%POS\t%INFO/ALLELEID\t%INFO/CLNHGVS\t%INFO/CLNREVSTAT\t%INFO/CLNSIG\n' -i 'INFO/CLNSIG="Likely_pathogenic" || INFO/CLNSIG="Pathogenic" || INFO/CLNSIG="Pathogenic/Likely_pathogenic"' -R non_chr_bed.tmp /reference_files/clinvar.vcf.gz > clinvar_annotation.txt
        rm non_chr_bed.tmp
    >>>
    output {
        File clinvar_annotation = "clinvar_annotation.txt"
    }
    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/annotatebed_hg38:v3"])
        memory: memory_gb + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
task concatenate_results {
    input{
        File clinvar_annotation
        File samtools_coverage_file
        File annotation_file
        String output_prefix
        Int memory_gb = 16
        Int disk_size = 32
        String? docker_override
    }
    command {
        python3 /scripts/aggregation_script.hg38.py --clinvar_file ~{clinvar_annotation} --samtools_coverage_file ~{samtools_coverage_file} --annotation_file ~{annotation_file} --output_prefix ~{output_prefix}
    }
    output {
        File integrated_annotation_file = "~{output_prefix}.integrated_annotation.txt"
    }
    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/annotatebed_hg38:v3"])
        memory: memory_gb + "GB"
        disks: "local-disk " + disk_size + " HDD"

    }
}

task GenerateGeneSummary {
    input {
       File CoverageFile
       File gene_bed
       File group_by_gene
       Int sample_fraction
       Int memory_gb = 32
       Int disk_size = 32
       Boolean generate_plot
    }

    command <<<

        Rscript /script/generate_gene_summary.R \
            --coverage_file ${CoverageFile} \
            --gene_bed ${gene_bed} \
            --grouped_by_gene ${group_by_gene} \
            --sample_fraction ${sample_fraction} \
            --plot ${if generate_plot then "true" else "false"}

    >>>

    output {
        File gene_base_count = "gene_base_count.txt"
        File? output_plot = "callable_fraction.png"
    }

    runtime {
        docker: "us.gcr.io/tag-public/samtools_r:v1"
        memory: memory_gb + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}