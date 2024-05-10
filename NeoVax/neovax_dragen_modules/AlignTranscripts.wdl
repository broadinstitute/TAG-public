version 1.0

import "./FormatUtil.wdl" as FormatUtil

struct Runtime {
  Int max_retries
  Int preemptible
  Int mem
  Int cpu
  Int boot_disk_size
  Int initial_disk_size
  String docker
}

workflow AlignTranscripts {
  input {
    String prefix

    # input format: BAM or FASTQ
    File? input_bam
    File? input_fastq1
    File? input_fastq2

    # Set 'true' to run the workflow in UMI aware mode
    Boolean umi_aware_mode = false

    # Set 'true' if UMIs need to be extracted from reads
    Boolean extract_umi = false
    # Read structure (e.g. 3M2S146T) to extract UMIs with fgbio
    # For more detials, see: https://fulcrumgenomics.github.io/fgbio/tools/latest/ExtractUmisFromBam.html
    String? read1_structure
    String? read2_structure

    # UMI-tools parameters
    String umi_tag_delimiter = "-"
    # UMI grouping method options: https://umi-tools.readthedocs.io/en/latest/reference/group.html
    String umi_grouping_method = "directional"
    String? umi_tools_extra_args

    # More STAR parameters in the STAR task
    File star_index
    Int star_cores = 8
    Int star_memory = 52

    # MarkDuplicates parameters
    Int max_records_in_ram = 500000
    Float sorting_collection_size_ratio = 0.25
    String tagging_policy = "DontTag"
    # Maximum offset between two duplicate clusters.
    # 100 (default) is appropriate for unpatterned,
    # 2500 recommended for patterned flowcells.
    Int optical_duplicate_pixel_distance = 100
    String? markduplicate_extra_args

    # workflow level runtime parameters
    Int preemptible = 2
    Int max_retries = 1
    Int additional_disk = 80
    Int boot_disk_size = 15
    Int mem = 16
    Int cpu = 2
  }

  Runtime standard_runtime = { "preemptible": preemptible,
                               "max_retries": max_retries,
                               "mem": mem,
                               "cpu": cpu,
                               "docker": "us.gcr.io/tag-public/neovax-tag-rnaseq:v1",
                               "boot_disk_size": boot_disk_size,
                               "initial_disk_size": additional_disk }
                               
  Runtime star_runtime = { "preemptible": preemptible,
                           "max_retries": max_retries,
                           "mem": star_memory,
                           "cpu": star_cores,
                           "docker": "us.gcr.io/tag-public/neovax-tag-rnaseq:v1",
                           "boot_disk_size": boot_disk_size,
                           "initial_disk_size": additional_disk*2 }

  if(defined(input_bam)) {
    call FormatUtil.SamToFastq as SamToFastq {
      input:
        prefix = prefix,
        input_bam = select_first([input_bam]),
        max_retries = standard_runtime.max_retries,
        preemptible = standard_runtime.preemptible,
        mem = standard_runtime.mem,
        cpu = standard_runtime.cpu,
        boot_disk_size = standard_runtime.boot_disk_size,
        disk_pad = standard_runtime.initial_disk_size
    }
  }

  File fastq1 = select_first([SamToFastq.fastq1, input_fastq1])
  File fastq2 = select_first([SamToFastq.fastq2, input_fastq2])
  call FormatUtil.FastqToSam as FastqToSam {
    input:
      prefix = prefix,
      fastq1 = fastq1,
      fastq2 = fastq2,
      mem = standard_runtime.mem,
      cpu = standard_runtime.cpu,
      boot_disk_size = standard_runtime.boot_disk_size,
      disk_pad = standard_runtime.initial_disk_size,
      preemptible = standard_runtime.preemptible,
      max_retries = standard_runtime.max_retries
  }

  if(extract_umi){
    call ExtractUMIs {
      input:
        prefix = prefix,
        unmapped_bam = FastqToSam.unmapped_bam,
        read1_structure = select_first([read1_structure]),
        read2_structure = select_first([read2_structure]),
        runtime_params = standard_runtime
    }
    call FormatUtil.SamToFastq as SamToFastqUMI {
      input:
        prefix = prefix,
        input_bam = ExtractUMIs.output_bam,
        max_retries = standard_runtime.max_retries,
        preemptible = standard_runtime.preemptible,
        mem = standard_runtime.mem,
        cpu = standard_runtime.cpu,
        boot_disk_size = standard_runtime.boot_disk_size,
        disk_pad = standard_runtime.initial_disk_size
    }
  }

  call FastQC as FirstReadsFastQC {
    input:
      fastq = select_first([SamToFastqUMI.fastq1, fastq1]),
      runtime_params = standard_runtime
  }
  call FastQC as SecondReadsFastQC {
    input:
      fastq = select_first([SamToFastqUMI.fastq2, fastq2]),
      runtime_params = standard_runtime
  }

  Boolean run_umi_aware_mode = if extract_umi then true else umi_aware_mode
  call STAR {
    input:
      prefix = prefix,
      input_bam = select_first([ExtractUMIs.output_bam, FastqToSam.unmapped_bam]),
      star_index = star_index,
      runtime_params = star_runtime,
      sort_transcriptome_bam = run_umi_aware_mode
  }

  if(run_umi_aware_mode) {
    call GroupByUMIs as GroupGenomeUMIs {
      input:
        prefix = prefix + ".Aligned.sortedByCoord.umiGrouped.out",
        input_bam = STAR.genome_bam,
        umi_tag_delimiter = umi_tag_delimiter,
        umi_grouping_method = umi_grouping_method,
        umi_tools_extra_args = umi_tools_extra_args,
        runtime_params = star_runtime
    }
    call GroupByUMIs as GroupTranscriptomeUMIs {
      input:
        prefix = prefix + ".Aligned.toTranscriptome.umiGrouped.out",
        input_bam = STAR.transcriptome_bam,
        umi_tag_delimiter = umi_tag_delimiter,
        umi_grouping_method = umi_grouping_method,
        umi_tools_extra_args = umi_tools_extra_args,
        runtime_params = star_runtime
    }
    call MarkDuplicates as MarkDuplicatesTranscriptome {
      input:
        prefix = prefix + ".Aligned.toTranscriptome.out",
        input_bam = GroupTranscriptomeUMIs.grouped_bam,
        is_transcriptome = true,
        remove_duplicates = true,
        run_umi_aware_mode = run_umi_aware_mode,
        max_records_in_ram = max_records_in_ram,
        sorting_collection_size_ratio = sorting_collection_size_ratio,
        tagging_policy = tagging_policy,
        optical_duplicate_pixel_distance = optical_duplicate_pixel_distance,
        markduplicate_extra_args = markduplicate_extra_args,
        runtime_params = star_runtime
    }
    call FormatTranscriptomeUMI {
      input:
        prefix = prefix + ".Aligned.toTranscriptome.out.md",
        input_bam = MarkDuplicatesTranscriptome.output_bam,
        runtime_params = star_runtime
    }
  }

  call MarkDuplicates as MarkDuplicatesGenome {
    input:
      prefix = prefix + ".Aligned.sortedByCoord.out.md",
      input_bam = select_first([GroupGenomeUMIs.grouped_bam, STAR.genome_bam]),
      is_transcriptome = false,
      remove_duplicates = false,
      run_umi_aware_mode = run_umi_aware_mode,
      max_records_in_ram = max_records_in_ram,
      sorting_collection_size_ratio = sorting_collection_size_ratio,
      tagging_policy = tagging_policy,
      optical_duplicate_pixel_distance = optical_duplicate_pixel_distance,
      markduplicate_extra_args = markduplicate_extra_args,
      runtime_params = star_runtime
  }
  
  output {
    File output_bam = MarkDuplicatesGenome.output_bam
    File output_bam_index = MarkDuplicatesGenome.output_bam_index
    File duplicate_metrics = MarkDuplicatesGenome.metrics
    
    File transcriptome_bam = select_first([FormatTranscriptomeUMI.output_bam, STAR.transcriptome_bam])
    File? transcriptome_duplicate_metrics = MarkDuplicatesTranscriptome.metrics

    File read1_fastqc_report = FirstReadsFastQC.fastqc_html
    File read1_fastqc_table = FirstReadsFastQC.fastqc_data
    File read2_fastqc_report = SecondReadsFastQC.fastqc_html
    File read2_fastqc_table = SecondReadsFastQC.fastqc_data    
  }
}

task ExtractUMIs {
  input {
    String prefix
    File unmapped_bam
    String read1_structure
    String read2_structure
    Runtime runtime_params
  }

  Int compute_mem = runtime_params.mem - 1
  Int disk_gb = ceil(size(unmapped_bam,"GB")*2.5) + runtime_params.initial_disk_size

  command {
    set -euo pipefail
    
    java "-Xmx${compute_mem}g" -jar /root/fgbio-1.5.0/fgbio.jar ExtractUmisFromBam \
      --input ~{unmapped_bam} \
      --output ~{prefix}.umi.u.bam \
      --molecular-index-tags "RX" \
      --read-structure ~{read1_structure} \
      --read-structure ~{read2_structure}
  }

  output {
    File output_bam = "~{prefix}.umi.u.bam"
  }

  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }

  meta {
    author: "Junko Tsuji"
  }
}

task FastQC {
  input {
    File fastq
    Runtime runtime_params
  }

  String fastq_basename = basename(fastq, ".fastq.gz")
  Int disk_gb = runtime_params.initial_disk_size + ceil(size(fastq,"GB"))

  command {
    set -euo pipefail

    fastqc ~{fastq} --threads ~{runtime_params.cpu} --extract --outdir .
    mv ~{fastq_basename}_fastqc/fastqc_data.txt ~{fastq_basename}_fastqc_data.txt
  }

  output {
    File fastqc_html = "~{fastq_basename}_fastqc.html"
    File fastqc_data = "~{fastq_basename}_fastqc_data.txt"
  }

  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries  
  }

  meta {
    author: "Junko Tsuji"
  }
}

task STAR {
  input {
    String prefix
    File input_bam
    File star_index
    Runtime runtime_params
    Boolean sort_transcriptome_bam

    # STAR related parameters
    Int outFilterMultimapNmax = 20
    Int alignSJoverhangMin = 8
    Int alignSJDBoverhangMin = 1
    Int outFilterMismatchNmax = 999
    Float outFilterMismatchNoverLmax = 0.1
    Int alignIntronMin = 20
    Int alignIntronMax = 1000000
    Int alignMatesGapMax = 1000000
    String outFilterType = "BySJout"
    Float outFilterScoreMinOverLread = 0.33
    Float outFilterMatchNminOverLread = 0.33
    Int limitSjdbInsertNsj = 1200000
    String outSAMstrandField = "intronMotif"
    String outFilterIntronMotifs = "None"
    String alignSoftClipAtReferenceEnds = "Yes"
    String quantMode = "TranscriptomeSAM GeneCounts"
    Int chimSegmentMin = 15
    Int chimJunctionOverhangMin = 15
    Int chimMainSegmentMultNmax = 1
    String chimOutType = "Junctions WithinBAM SoftClip"
    String outSAMattributes = "NH HI AS nM NM ch"

    # SJ.out.tab file from 1st pass. Only one pass will be run with this option
    File? sjdbFileChrStartEnd
    
    # if fragment size is short, and there are many overlapping reads, increase this to 20-30
    String alignEndsProtrude = "10 ConcordantPair"
    
    String? star_extra_args
    
    # Variant VCF file with germline variants. Only SNPs are supported at the moment,
    # Each variant is expected to have a genotype with two alleles - e.g. 0/1
    File? varVCFfile
  }

  Int disk_gb = runtime_params.initial_disk_size + ceil(size(varVCFfile,"GB") * 2.5) +
                ceil(size(input_bam,"GB") + size(star_index,"GB") * 4.5)
  
  command {
    set -euo pipefail

    # extract index
    mkdir star_index
    tar -xvvf ~{star_index} -C star_index --strip-components=1

    # splice junction parameters
    sjdb_params=" --twopassMode Basic "
    if [[ -n "~{sjdbFileChrStartEnd}" ]] ; then
      sjdb_params=" --sjdbFileChrStartEnd ~{sjdbFileChrStartEnd} "
    fi

    # chimeric read parameters
    chimera_params=" --chimSegmentMin ~{chimSegmentMin} "
    if [[ ~{chimSegmentMin} -gt 0 ]] ; then
      chimera_params=" --chimJunctionOverhangMin ~{chimJunctionOverhangMin} "$chimera_params
      chimera_params=" --chimMainSegmentMultNmax ~{chimMainSegmentMultNmax} "$chimera_params
      chimera_params=" --chimOutJunctionFormat 0 --chimOutType ~{chimOutType} "$chimera_params
    fi

    # variant detection & WASP mode
    sam_attribute_params=" --outSAMattributes ~{outSAMattributes} "
    var_params="~{varVCFfile}"
    if [[ -n "~{varVCFfile}" ]] ; then
      if [[ `basename ~{varVCFfile} .gz` != `basename ~{varVCFfile}` ]] ; then
         var_params=`basename ~{varVCFfile} .gz`
         mv ~{varVCFfile} . && gunzip `basename ~{varVCFfile}`
      fi
      var_params=" --varVCFfile "$varparams" --waspOutputMode SAMtag"
      sam_attribute_params=$sam_attribute_params" vA vG vW"
    fi

    STAR --runMode alignReads \
      --runThreadN ~{runtime_params.cpu} \
      --genomeDir star_index \
      --readFilesIn ~{input_bam} --readFilesType SAM PE \
      --readFilesCommand samtools view -h \
      --outFileNamePrefix "~{prefix}." \
      ~{"--outFilterMultimapNmax " + outFilterMultimapNmax} \
      ~{"--alignSJoverhangMin " + alignSJoverhangMin} \
      ~{"--alignSJDBoverhangMin " + alignSJDBoverhangMin} \
      ~{"--outFilterMismatchNmax " + outFilterMismatchNmax} \
      ~{"--outFilterMismatchNoverLmax " + outFilterMismatchNoverLmax} \
      ~{"--alignIntronMin " + alignIntronMin} \
      ~{"--alignIntronMax " + alignIntronMax} \
      ~{"--alignMatesGapMax " + alignMatesGapMax} \
      ~{"--outFilterType " + outFilterType} \
      ~{"--outFilterScoreMinOverLread " + outFilterScoreMinOverLread} \
      ~{"--outFilterMatchNminOverLread " + outFilterMatchNminOverLread} \
      ~{"--limitSjdbInsertNsj " + limitSjdbInsertNsj} \
      ~{"--outSAMstrandField " + outSAMstrandField} \
      ~{"--outFilterIntronMotifs " + outFilterIntronMotifs} \
      ~{"--alignSoftClipAtReferenceEnds " + alignSoftClipAtReferenceEnds} \
      ~{"--alignEndsProtrude " + alignEndsProtrude} \
      ~{"--quantMode " + quantMode} \
      --outSAMtype BAM Unsorted \
      --outSAMunmapped Within \
      --genomeLoad NoSharedMemory \
      --outSAMattrRGline "ID:RG_~{prefix}" "SM:~{prefix}" "LB:~{prefix}" "PL:ILLUMINA" "PU:barcode_~{prefix}" \
      $sjdb_params \
      $sam_attribute_params \
      $chimera_params \
      $var_params

    # sort genome aligned BAM
    samtools sort --threads ~{runtime_params.cpu} -o ~{prefix}.Aligned.sortedByCoord.out.bam ~{prefix}.Aligned.out.bam
    rm ~{prefix}.Aligned.out.bam
    samtools index ~{prefix}.Aligned.sortedByCoord.out.bam

    # sort transcritome BAM if UMIs are attached
    if [[ "~{sort_transcriptome_bam}" = "true" ]] && [[ -f "~{prefix}.Aligned.toTranscriptome.out.bam " ]] ; then
      samtools sort --threads ~{runtime_params.cpu} -o ~{prefix}.tr.bam ~{prefix}.Aligned.toTranscriptome.out.bam
      mv ~{prefix}.tr.bam ~{prefix}.Aligned.toTranscriptome.out.bam
    fi

    # sort chimeric read BAM
    if [[ -f "~{prefix}.Chimeric.out.sam" ]] ; then
      samtools sort --threads ~{runtime_params.cpu} -o ~{prefix}.Chimeric.out.sorted.bam ~{prefix}.Chimeric.out.sam
      samtools index ~{prefix}.Chimeric.out.sorted.bam
    fi

    # compress splice and chimeric junction files
    gzip ~{prefix}.SJ.out.tab
    mv ~{prefix}._STARpass1/SJ.out.tab ~{prefix}.SJ.pass1.out.tab && gzip ~{prefix}.SJ.pass1.out.tab
    
    if [[ -f "~{prefix}.Chimeric.out.junction" ]] ; then
      gzip ~{prefix}.Chimeric.out.junction
    fi

    # compress read count file
    if [[ -f "~{prefix}.ReadsPerGene.out.tab" ]] ; then
      gzip ~{prefix}.ReadsPerGene.out.tab
    fi

    # placeholders for optional outputs
    touch ~{prefix}.Aligned.toTranscriptome.out.bam
    touch ~{prefix}.Chimeric.out.sorted.bam
    touch ~{prefix}.Chimeric.out.sorted.bam.bai
    touch ~{prefix}.ReadsPerGene.out.tab.gz
  }

  output {
    File genome_bam = "~{prefix}.Aligned.sortedByCoord.out.bam"
    File genome_bam_index = "~{prefix}.Aligned.sortedByCoord.out.bam.bai"
    File transcriptome_bam = "~{prefix}.Aligned.toTranscriptome.out.bam"
    File junctions = "~{prefix}.SJ.out.tab.gz"
    File junctions_pass1 = "~{prefix}.SJ.pass1.out.tab.gz"
    File chimeric_junctions = "~{prefix}.Chimeric.out.junction.gz"
    File chimeric_bam = "~{prefix}.Chimeric.out.sorted.bam"
    File chimeric_bam_index = "~{prefix}.Chimeric.out.sorted.bam.bai"
    File read_counts = "~{prefix}.ReadsPerGene.out.tab.gz"
  }

  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }

  meta {
    author: "Junko Tsuji"
  }
}

task GroupByUMIs {
  input {
    String prefix
    File input_bam
    String umi_tag_delimiter
    String umi_grouping_method
    String? umi_tools_extra_args    
    Runtime runtime_params
  }

  Int disk_gb = runtime_params.initial_disk_size + ceil(size(input_bam,"GB") * 4.5)

  command {
    set -euo pipefail

    # Group UMIs with UMI tools
    umi_tools group -I ~{input_bam} --paired --no-sort-output \
      --output-bam --stdout ~{prefix}.grouped.bam \
      --extract-umi-method tag --umi-tag "RX" --unmapped-reads use \
      --umi-tag-delimiter ~{umi_tag_delimiter} --umi-group-tag "BX" \
      --method ~{umi_grouping_method} ~{umi_tools_extra_args}

    # Sort generated BAM in the query name order
    samtools sort --threads ~{runtime_params.cpu} -n \
      -o ~{prefix}.bam ~{prefix}.grouped.bam
  }

  output {
    File grouped_bam = "~{prefix}.bam"
  }

  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }

  meta {
    author: "Junko Tsuji"
  }
}


task MarkDuplicates {
  input {
    String prefix
    File input_bam
    Boolean is_transcriptome
    Boolean run_umi_aware_mode
    Boolean remove_duplicates

    Int max_records_in_ram
    Float sorting_collection_size_ratio
    String tagging_policy
    Int optical_duplicate_pixel_distance
    String? markduplicate_extra_args

    Runtime runtime_params    
  }

  String input_bam_sort_order = if run_umi_aware_mode then "queryname" else "coordinate"
  Int disk_gb = runtime_params.initial_disk_size + ceil(size(input_bam,"GB")*2.5)
  Int compute_mem = runtime_params.mem - 1
  
  command {
    java "-Xmx~{compute_mem}g" -jar /root/picard-2.26.10/picard.jar MarkDuplicates \
      -I ~{input_bam} \
      -O ~{prefix}.bam \
      -M ~{prefix}.duplicate_metrics.txt \
      --ASSUME_SORT_ORDER ~{input_bam_sort_order} \
      ~{true="--READ_ONE_BARCODE_TAG BX" false="" run_umi_aware_mode} \
      ~{true="--REMOVE_DUPLICATES" false="" remove_duplicates} \
      --MAX_RECORDS_IN_RAM ~{max_records_in_ram} \
      --SORTING_COLLECTION_SIZE_RATIO ~{sorting_collection_size_ratio} \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE ~{optical_duplicate_pixel_distance} \
      --TAGGING_POLICY ~{tagging_policy} \
      ~{markduplicate_extra_args}

    # Only sort query-name-sorted 'genome' BAM by coordinates
    if [[ "~{input_bam_sort_order}" = "queryname" ]] && [[ "~{is_transcriptome}" = "false" ]] ; then
      samtools sort --threads ~{runtime_params.cpu} -o ~{prefix}.sorted.bam ~{prefix}.bam
      mv ~{prefix}.sorted.bam ~{prefix}.bam
    fi
    # Create BAM index for coordinate sorted BAM
    if [[ "~{is_transcriptome}" = "false" ]] ; then
      samtools index ~{prefix}.bam
    fi
    # touch the index file to avoid error
    touch ~{prefix}.bam.bai
  }

  output {
    File output_bam = "~{prefix}.bam"
    File output_bam_index = "~{prefix}.bam.bai"
    File metrics = "~{prefix}.duplicate_metrics.txt"
  }

  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }

  meta {
    author: "Junko Tsuji"
  }
}

task FormatTranscriptomeUMI {
  input {
    String prefix
    File input_bam
    Runtime runtime_params
  }

  Int disk_gb = runtime_params.initial_disk_size + ceil(size(input_bam,"GB")*2.5)

  command {
    umi_tools prepare-for-rsem --tags UG,BX,RX \
      -I ~{input_bam} --stdout ~{prefix}.bam
  }
  
  output {
    File output_bam = "~{prefix}.bam"
  }
  
  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }

  meta {
    author: "Junko Tsuji"
  }
}