# This WDL takes a sample set and checks the coverage and the number
# of variants detected in specific genes. The gene BED file needs to
# have gene names in the name column (the 4th column).

workflow CollectGeneStatistics {
   Array[File] input_bams
   Array[File] input_bams_index
   Array[File] variant_vcfs
   Array[File] variant_vcfs_index
   String sample_set_name

   String funcotator_docker
   File funcotator_annotation_tar_gz
   
   File gene_bed
   File reference
   File reference_index
   File reference_dict
   Int ref_size = ceil(size(reference,"GB") + size(reference_index,"GB") + size(reference_dict,"GB"))
   
   Int? disk_pad = 20

   call BedToIntervalList {
      input:
         gene_bed = gene_bed,
         reference_dict = reference_dict,
         disk_pad = disk_pad
   }

   scatter (input_bam in zip(input_bams, input_bams_index)) {
      Int bam_size = ceil(size(input_bam.left,"GB") + size(input_bam.right,"GB"))
      call CollectGeneCoverage {
         input: 
            interval_list = BedToIntervalList.interval_list,
            gene_bed = gene_bed,
            reference = reference,
            reference_index = reference_index,
            reference_dict = reference_dict,
            input_bam = input_bam.left,
            input_bam_index = input_bam.right,
            out_basename = basename(input_bam.left, ".bam"),
            disk_size = bam_size * 2 + ref_size + disk_pad
      }
   }
   
   scatter (input_vcf in zip(variant_vcfs, variant_vcfs_index)) {
      Int vcf_size = ceil(size(input_vcf.left,"GB") + size(input_vcf.right,"GB"))
      Boolean is_not_compressed = basename(input_vcf.left, ".vcf.gz") == basename(input_vcf.left)
      String out_basename = if is_not_compressed then basename(input_vcf.left, ".vcf") else basename(input_vcf.left, ".vcf.gz")

      call Funcotator {
         input:
            reference = reference,
            reference_index = reference_index,
            reference_dict = reference_dict,
            input_vcf = input_vcf.left,
            input_vcf_index = input_vcf.right,
            out_basename = out_basename,
            reference_version = "hg19",
            data_sources_tar_gz = funcotator_annotation_tar_gz,
            tumor_name = out_basename,
            funcotator_docker = funcotator_docker,
            disk_size = ceil(size(funcotator_annotation_tar_gz, "GB") * 8) + ref_size + vcf_size + disk_pad
      }

      call CollectGeneMutations {
         input:
            gene_bed = gene_bed,
            input_vcf = input_vcf.left,
            input_vcf_index = input_vcf.right,
            input_maf = Funcotator.funcotated_maf,
            out_basename = if is_not_compressed then basename(input_vcf.left, ".vcf") else basename(input_vcf.left, ".vcf.gz"),
            disk_size = vcf_size + ceil(size(Funcotator.funcotated_maf,"GB")) + disk_pad
      }
   }

   call GatherStatistics {
      input:
         out_basename = sample_set_name,
         coverage_fractions = CollectGeneCoverage.coverage_fraction,
         variant_stats = CollectGeneMutations.variant_stats
   }

   output {
      File variant_stats_set = GatherStatistics.variant_stats_set
      File coverage_fraction_set = GatherStatistics.coverage_fraction_set
   }
}


task BedToIntervalList {
   File gene_bed
   File reference_dict
   String out_basename = basename(gene_bed, ".bed")
   Int? disk_pad
   
   Int disk_size = ceil(size(gene_bed, "GB")) * 2 + disk_pad

   command {
      set -e
      
      /gatk/gatk BedToIntervalList \
        -SD ${reference_dict} \
        -I ${gene_bed} \
        -O "${out_basename}.interval_list"
   }
   runtime {
      docker: "us.gcr.io/broad-gatk/gatk:4.1.9.0"
      memory: "4GB"
      disks: "local-disk " + disk_size + " HDD"
      maxRetries: 1
      preemptible: 3
   }
   output {
      File interval_list = "${out_basename}.interval_list"
   }
}

task CollectGeneCoverage {
   File interval_list
   File gene_bed
   File reference
   File reference_index
   File reference_dict
   
   File input_bam
   File input_bam_index
   String out_basename

   Int? min_mapping_quality = 20
   Int? min_base_quality = 20
   String? extra_arguments = "--countType COUNT_FRAGMENTS_REQUIRE_SAME_BASE"
   Array[Int]? coverage_cutoffs = [100, 500, 1000]

   Int disk_size
   Int? preemptible_attempts = 3
   Int? mem = 8

   command <<<
       set -e

       # compute per-base coverage
       java -jar /usr/tag/GATK36.jar -T DepthOfCoverage \
         -L ${interval_list} \
         -I ${input_bam} \
         -R ${reference} \
         -o ${out_basename}.depth \
         --omitPerSampleStats \
         -mmq ${min_mapping_quality} \
         -mbq ${min_base_quality} \
         ${extra_arguments}

       # covert per-base coverage into bed and link the positions to gene names
       grep -v ^Locus "${out_basename}.depth" | \
       awk '{split($1,pos,":"); print pos[1],pos[2]-1,pos[2],$2}' OFS="\t" | \
       bedtools intersect -a "${gene_bed}" -b stdin -wo | \
       cut -f1-4,8 > "${out_basename}_cov.bed"

       # compute total bases of genes of interest
       awk '{arr[$4]+=$3-$2} END{for(i in arr){print i"\t"arr[i]}}' "${gene_bed}" | sort -k1,1 > total_gene_bases.out
       
       # calculate covered base fraction at a specific cutoff
       for CUTOFF in ${sep=" " coverage_cutoffs};
       do
        awk -v cov=$CUTOFF \
            -v gene_str=`cut -f1 total_gene_bases.out | tr "\n" ","` \
            ' BEGIN{   split(gene_str, genes, ",");
                       for(i in genes){ if(length(genes[i]) > 0){ arr[genes[i]] = 0 } }
                   }
             {  if($5 >= cov){ arr[$4]+=1 } }
              END{ for(i in arr){ print i"\t"arr[i] } }' "${out_basename}_cov.bed" | sort -k1,1 > sample_gene_bases_$CUTOFF.out
        
        # calculate covered base fraction of genes
        paste sample_gene_bases_$CUTOFF.out total_gene_bases.out | \
        awk -v cov=$CUTOFF -v sample=${out_basename} \
           '{print sample,$1,$2/$4,"cutoff="cov"X"}' OFS="\t" >> "${out_basename}.coverage_fraction.txt"
       done
   >>>
   runtime {
      docker: "us.gcr.io/tag-team-160914/tag-tools:1.0.0"
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + "GB"
      preemptible: preemptible_attempts
      maxRetries: 1      
   }
   output {
      File coverage_fraction = "${out_basename}.coverage_fraction.txt"
   }
}

task CollectGeneMutations {
   File gene_bed
   File input_vcf
   File input_vcf_index
   File input_maf
   String out_basename

   Int disk_size
   Int? preemptible_attempts = 3
   Int? mem = 4
   
   command <<<
      # python code snippet to create BED from MAF
      ## ==== python code block starts ====
PYTHON_SRC="
columns_to_keep = ['Chromosome', 'Start_Position', 'End_Position',
                   'Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification',
                   'cDNA_Change', 'Protein_Change', 'tumor_f', 't_alt_count', 't_ref_count']
fo = open('${out_basename}_maf.bed', 'w')
for line in open('${input_maf}'):
    if line.startswith('#'):
        continue
    line = line.rstrip('\n').split('\t')    
    if line[0] == 'Hugo_Symbol':
        index = [line.index(key) for key in columns_to_keep]
        continue
    line = [line[i] if line[i] else 'NA' for i in index]
    line[1] = str(int(line[1])-1)  # adjusting start offset 
    fo.write('\t'.join(line) + '\n')
fo.close()
"
      ## ==== python code block ends ====

      # convert MAF to BED
      python -c "$PYTHON_SRC"

      # intersect variants against gene intervals
      bedtools intersect -a ${out_basename}_maf.bed -b ${gene_bed} | \
      bedtools intersect -a stdin -b ${input_vcf} -wo | \
      awk 'BEGIN{ print "sample\tgene\tposition\tfilter_flag\tvariant_class\tcDNA_change\tprotein_change\ttumor_f\tt_alt_count\tt_ref_count" }             
           {  print $4,$5,$1":"$2+1"-"$3,$18,$6,$7,$8,$9,$10,$11 }' OFS="\t" > "${out_basename}.variant_stats.txt"
   >>>
   runtime {
      docker: "us.gcr.io/tag-team-160914/tag-tools:1.0.0"
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + "GB"
      preemptible: preemptible_attempts
      maxRetries: 1
   }
   output {
      File variant_stats = "${out_basename}.variant_stats.txt"
   }
}


task GatherStatistics {
   String out_basename
   Array[File] variant_stats
   Array[File] coverage_fractions

   command <<<
      set -e
      cat ${sep=" " variant_stats} | grep -v ^sample | \
      cat <(grep ^sample ${variant_stats[0]}) - >> ${out_basename}.variant_stats.txt
      cat ${sep=" " coverage_fractions} >> ${out_basename}.coverage_fraction.txt
   >>>
   runtime {
      docker: "us.gcr.io/tag-team-160914/tag-tools:1.0.0"
      disks: "local-disk 20 HDD"
      memory: "4GB"
      maxRetries: 1      
      preemptible: 3
   }
   output {
       File variant_stats_set = "${out_basename}.variant_stats.txt"
       File coverage_fraction_set = "${out_basename}.coverage_fraction.txt"       
   }
}


task Funcotator {
   File reference
   File reference_index
   File reference_dict
   File input_vcf
   File input_vcf_index
   String out_basename
   String reference_version
   File data_sources_tar_gz
   String tumor_name
   String? normal_name
   String? sequencing_center
   String? sequence_source   
   String? extra_args

   String funcotator_docker
   String? gatk_override
   Int? mem = 4
   Int command_mem = mem * 1000 - 500
   Int disk_size

   command <<<
      set -e
      export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

      # Extract annotation database
      mkdir datasources_dir
      tar zxvf ${data_sources_tar_gz} -C datasources_dir --strip-components 1
      DATA_SOURCES_FOLDER="$PWD/datasources_dir"

      # Run Funcotator
      gatk --java-options "-Xmx${command_mem}m" Funcotator \
           --data-sources-path $DATA_SOURCES_FOLDER \
           --ref-version ${reference_version} \
           --output-file-format MAF \
           -R ${reference} \
           -V ${input_vcf} \
           -O "${out_basename}.funcotated.maf" \
           --annotation-default tumor_barcode:${default="Unknown" tumor_name} \
           --annotation-default normal_barcode:${default="Unknown" normal_name} \
           --annotation-default Center:${default="Unknown" sequencing_center} \
           --annotation-default source:${default="Unknown" sequence_source} \
           ${extra_args}
     >>>

    runtime {
        docker: funcotator_docker
        bootDiskSizeGb: 12
        memory: mem + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
        maxRetries: 1
    }

     output {
         File funcotated_maf = "${out_basename}.funcotated.maf"
     }
}
