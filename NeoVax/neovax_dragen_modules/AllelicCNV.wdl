version 1.0

struct Runtime {
  Int max_retries
  Int preemptible
  Int mem
  Int cpu
  Int boot_disk_size
  Int initial_disk_size
  String docker
}

workflow AllelicCNV {
  input {
    String pair_name
    
    File tumor_bam
    File tumor_bam_index

    # this version allows to do tumor-only CNV calls
    File? normal_bam
    File? normal_bam_index

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File cnv_pon
    File common_snp_list

    # workflow level runtime parameters
    Int preemptible = 2
    Int max_retries = 1
    Int additional_disk = 50
    Int boot_disk_size = 15
    Int mem = 16
    Int cpu = 1
  }

  Runtime standard_runtime = { "preemptible": preemptible,
                               "max_retries": max_retries,
                               "mem": mem,
                               "cpu": cpu,
                               "docker": "us.gcr.io/tag-public/neovax-tag-cga:v1",
                               "boot_disk_size": boot_disk_size,
                               "initial_disk_size": additional_disk }
    
  call GatkACNV {
    input:
      pair_name = pair_name,
      tumor_bam = tumor_bam,
      tumor_bam_index = tumor_bam_index,
      normal_bam = normal_bam,
      normal_bam_index = normal_bam_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      cnv_pon = cnv_pon,
      common_snp_list = common_snp_list,
      runtime_params = standard_runtime
  }

  call AllelicCapSeg {
    input:
      pair_name = pair_name,
      input_seg = GatkACNV.output_seg,
      post_tn_coverage = GatkACNV.post_tn_coverage,
      tumor_hets = GatkACNV.tumor_hets,
      runtime_params = standard_runtime
  }

  output {
    File output_seg = GatkACNV.output_seg
    File raw_coverage = GatkACNV.raw_coverage
    File pre_tn_coverage = GatkACNV.pre_tn_coverage
    File post_tn_coverage = GatkACNV.post_tn_coverage
    File normal_hets = GatkACNV.normal_hets
    File tumor_hets = GatkACNV.tumor_hets
    File copy_ratio_plot = GatkACNV.copy_ratio_plot
    File copy_ratio_plot_lim4 = GatkACNV.copy_ratio_plot_lim4
    Float pre_tn_mad = GatkACNV.pre_tn_mad
    Float post_tn_mad = GatkACNV.post_tn_mad
	File alleliccapseg_plot = AllelicCapSeg.alleliccapseg_plot
	File alleliccapseg_tsv = AllelicCapSeg.alleliccapseg_tsv
	Float alleliccapseg_skew = AllelicCapSeg.alleliccapseg_skew
  }
}


task GatkACNV {
  input {
    String pair_name
    File tumor_bam
    File tumor_bam_index
    File? normal_bam
    File? normal_bam_index

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File cnv_pon
    File common_snp_list

    Runtime runtime_params
  }

  Int tumor_size = ceil(size(tumor_bam,"GB") + size(tumor_bam_index,"GB"))
  Int normal_size = ceil(size(normal_bam,"GB") + size(normal_bam_index,"GB"))
  Int ref_size = ceil(size(ref_fasta,"GB") + size(ref_fasta_index,"GB") + size(ref_dict,"GB"))
  
  Int disk_gb = tumor_size + normal_size + ref_size + runtime_params.initial_disk_size
  Int compute_mem = runtime_params.mem - 1

  command <<<
    set -euxo pipefail
         
    EFFECTIVE_TARGETS="padded_bed.bed"
    python /usr/local/bin/get_targets.py ~{cnv_pon} $EFFECTIVE_TARGETS

    # Coverage collection
    /usr/local/jre1.8.0_73/bin/java "-Xmx~{compute_mem}g" -Djava.library.path=/usr/lib/jni/ \
      -jar /root/gatk-protected.jar CalculateTargetCoverage \
      --output ~{pair_name}.coverage.tsv \
      --groupBy SAMPLE \
      --transform PCOV \
      --targets $EFFECTIVE_TARGETS \
      --targetInformationColumns FULL \
      --input ~{tumor_bam} \
      --reference ~{ref_fasta} \
      --interval_set_rule UNION \
      --interval_padding 0 \
      --secondsBetweenProgressUpdates 10.0 \
      --disableSequenceDictionaryValidation false \
      --createOutputBamIndex true \
      --help false \
      --version false \
      --verbosity INFO \
      --QUIET false

    # Normalization
    /usr/local/jre1.8.0_73/bin/java "-Xmx~{compute_mem}g" -Djava.library.path=/usr/lib/jni/ \
      -jar /root/gatk-protected.jar NormalizeSomaticReadCounts \
      --input ~{pair_name}.coverage.tsv \
      --targets $EFFECTIVE_TARGETS \
      --panelOfNormals ~{cnv_pon} \
      --factorNormalizedOutput ~{pair_name}.fnt.tsv \
      --tangentNormalized ~{pair_name}.tn.tsv \
      --betaHatsOutput ~{pair_name}.betaHats.tsv \
      --preTangentNormalized ~{pair_name}.pre_tn.tsv \
      --help false \
      --version false \
      --verbosity INFO \
      --QUIET false
      
    # Segmentation by tangent-normalized coverage
    /usr/local/jre1.8.0_73/bin/java "-Xmx~{compute_mem}g" -Djava.library.path=/usr/lib/jni/ \
      -jar /root/gatk-protected.jar PerformSegmentation \
      --tangentNormalized ~{pair_name}.tn.tsv \
      --output ~{pair_name}.seg \
      --log2Input true \
      --alpha 0.01 \
      --nperm 10000 \
      --pmethod HYBRID \
      --minWidth 2 \
      --kmax 25 \
      --nmin 200 \
      --eta 0.05 \
      --trim 0.025 \
      --undoSplits NONE \
      --undoPrune 0.05 \
      --undoSD 3 \
      --help false \
      --version false \
      --verbosity INFO \
      --QUIET false

    /usr/local/jre1.8.0_73/bin/java "-Xmx~{compute_mem}g" -Djava.library.path=/usr/lib/jni/ \
      -jar /root/gatk-protected.jar CallSegments \
      --tangentNormalized ~{pair_name}.tn.tsv \
      --segments ~{pair_name}.seg \
      --output ~{pair_name}.called \
      --legacy false \
      --help false \
      --version false \
      --verbosity INFO \
      --QUIET false

    # Run GetHetCoverage for tumor-normal analysis
    # Run GetBayesianHetCoverage for tumor-only analysis
    if [[ -f "~{normal_bam}" ]]; then
      /usr/local/jre1.8.0_73/bin/java "-Xmx~{compute_mem}g" -Djava.library.path=/usr/lib/jni/ \
        -jar /root/gatk-protected.jar GetHetCoverage \
        --reference ~{ref_fasta} \
        --normal ~{normal_bam} \
        --tumor ~{tumor_bam} \
        --snpIntervals ~{common_snp_list} \
        --normalHets ~{pair_name}.normal.hets.tsv \
        --tumorHets ~{pair_name}.tumor.hets.tsv \
        --pvalueThreshold 0.05 \
        --help false \
        --version false \
        --verbosity INFO \
        --QUIET false \
        --VALIDATION_STRINGENCY LENIENT
    else
      /usr/local/jre1.8.0_73/bin/java "-Xmx~{compute_mem}g" -Djava.library.path=/usr/lib/jni/ \
        -jar /root/gatk-protected.jar GetBayesianHetCoverage \
        --reference ~{ref_fasta} \
        --tumor ~{tumor_bam} \
        --snpIntervals ~{common_snp_list} \
        --tumorHets ~{pair_name}.tumor.hets.tsv \
        --hetCallingStringency 30 \
        --readDepthThreshold 15 \
        --help false \
        --version false \
        --verbosity INFO \
        --QUIET false \
        --VALIDATION_STRINGENCY LENIENT
        
        # create a null normal hets tsv file to avoid an error
        echo null > ~{pair_name}.normal.hets.tsv
    fi
        
    mkdir -v plotting
    echo "plotting segmented copy ratio"

    /usr/local/jre1.8.0_73/bin/java "-Xmx~{compute_mem}g" -Djava.library.path=/usr/lib/jni/ \
      -jar /root/gatk-protected.jar PlotSegmentedCopyRatio \
      --tangentNormalized ~{pair_name}.tn.tsv \
      --preTangentNormalized ~{pair_name}.pre_tn.tsv \
      --sequenceDictionaryFile ~{ref_dict} \
      --outputPrefix ~{pair_name} \
      --segments ~{pair_name}.seg \
      --output plotting/ \
      --log2Input true \
      --help false \
      --version false \
      --verbosity INFO \
      --QUIET false

    # Extract MAD values
    # MAD before tangent-normalization
    grep -v ^"#" ~{pair_name}.pre_tn.tsv | grep -v ^contig | cut -f5 | \
    R --slave -e 'median(abs(diff(2**scan("stdin"))))' | cut -f2 -d" " > pre_tn_mad.txt
    # MAD after tangen-normalization
    grep -v ^"#" ~{pair_name}.tn.tsv | grep -v ^contig | cut -f5 | \
    R --slave -e 'median(abs(diff(2**scan("stdin"))))' | cut -f2 -d" " > post_tn_mad.txt
  >>>

  output {
    File output_seg = "~{pair_name}.seg"
    File raw_coverage = "~{pair_name}.coverage.tsv"
    File pre_tn_coverage = "~{pair_name}.pre_tn.tsv"
    File post_tn_coverage = "~{pair_name}.tn.tsv"
    File normal_hets = "~{pair_name}.normal.hets.tsv"
    File tumor_hets = "~{pair_name}.tumor.hets.tsv"
    File copy_ratio_plot = "plotting/~{pair_name}_Before_After.png"
    File copy_ratio_plot_lim4 = "plotting/~{pair_name}_Before_After_CR_Lim_4.png"
    Float pre_tn_mad = read_float("pre_tn_mad.txt")
    Float post_tn_mad = read_float("post_tn_mad.txt")
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

  # modified the 'gatk_acnv_only' task in the CGA pipeline
  meta {
    author: "Junko Tsuji"
  } 
}

task AllelicCapSeg {
  input {
    String pair_name
    File post_tn_coverage
    File tumor_hets
    File input_seg
    Runtime runtime_params
  }

  Int disk_gb = ceil(size(post_tn_coverage,"GB") + size(tumor_hets,"GB") +
                     size(input_seg,"GB")) + runtime_params.initial_disk_size

  command <<<
    # Remove header and reorder columns
    grep -v ^"#" ~{post_tn_coverage} | \
    awk -F "\t" '{print $4,$1,$2,$3,$NF}' OFS="\t" > post_tn_no_header.tsv

    # Reformat tumor het file and insert dummy bases
    HD="Chromosome\tStart_position\tReference_Allele\tTumor_Seq_Allele1\ti_t_ref_count\ti_t_alt_count"
    echo -e $HD > tumor_hets_reformatted.tsv
    grep -v ^CONTIG ~{tumor_hets} | \
    awk -F "\t" '{print $1,$2,"A","T",$3,$4}' OFS="\t" >> tumor_hets_reformatted.tsv

    # Re-scale segment mean in seg file
    python <<CODE
    import math
    fo = open('seg_rescaled.tsv', 'w')
    for line in open('~{input_seg}'):
        if line.startswith('Sample'):
            fo.write(line)
            continue
        line = line.strip('\n').split('\t')
        segment_mean = math.log(float(line[-1]))/math.log(2)
        fo.write('\t'.join(line[:-1]+[str(segment_mean)])+'\n')
    fo.close()
    CODE

    Rscript /opt/Allelic_CapSeg/AllelicCapseg_cli.R \
      --SID=~{pair_name} \
      --capseg.probe.fn=post_tn_no_header.tsv \
      --capseg.seg.fn=seg_rescaled.tsv \
      --germline.het.fn=tumor_hets_reformatted.tsv \
      --drop.x=FALSE \
      --drop.y=TRUE \
      --seg.merge.thresh=0.5 \
      --min.seg.size=3 \
      --verbose=TRUE \
      --base.output.dir=. \
      --initial.merge=TRUE \
      --split.merge=TRUE \
      --outlier.thresh=0.005 \
      --working.dir=/opt/Allelic_CapSeg/
  >>>
  
  output {
    File alleliccapseg_tsv = "results/~{pair_name}.tsv"
    File alleliccapseg_plot = "plots/~{pair_name}/~{pair_name}_FullGenome.png"
    Float alleliccapseg_skew = read_float("results/~{pair_name}.skew")
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

  # modified the 'gatk_acnv_only' task in the CGA pipeline
  meta {
    author: "Junko Tsuji"
  } 
}