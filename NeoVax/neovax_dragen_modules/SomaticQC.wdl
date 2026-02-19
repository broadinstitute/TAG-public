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

workflow SomaticQC {
  input {
    String pair_name
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File target_intervals

    # workflow tries to find SNVs and INDELs from callstats and VCF files
    # if callstats file is supplied
    File? input_callstats
    File input_vcf
    File input_seg
    File tumor_hets
    File normal_hets

    # deTiN parameters
    File exac_pickle
    Float TiN_prior = 0.2
    Float mutation_prior = 0.05

    # ContEst parameters
    File snp6_bed
    File hapmap_vcf

    # Orientation bias Q-value
    File pre_adapter_detail_metrics

    # workflow level runtime parameters
    Int preemptible = 2
    Int max_retries = 1
    Int additional_disk = 50
    Int boot_disk_size = 15
    Int mem = 8
    Int cpu = 1
  }

  Runtime standard_runtime = { "preemptible": preemptible,
                                "max_retries": max_retries,
                                "mem": mem,
                                "cpu": cpu,
                                "docker": "us.gcr.io/tag-public/neovax-tag-tumorqc:v1",
                                "boot_disk_size": boot_disk_size,
                                "initial_disk_size": additional_disk }

  Int ref_size = ceil(size(ref_fasta,"GB") + size(ref_fasta_index,"GB") + size(ref_dict,"GB"))
  Int tumor_size = ceil(size(tumor_bam,"GB") + size(tumor_bam_index,"GB"))
  Int normal_size = ceil(size(normal_bam,"GB") + size(normal_bam_index,"GB"))

  call FormatUtil.GetSampleID as GetTumorID {
    input:
      bam = tumor_bam,
      bam_index = tumor_bam_index,
      max_retries = standard_runtime.max_retries,
      preemptible = standard_runtime.preemptible,
      mem = standard_runtime.mem,
      cpu = standard_runtime.cpu,
      boot_disk_size = standard_runtime.boot_disk_size,
      disk_pad = standard_runtime.initial_disk_size
  }

  call FormatUtil.GetSampleID as GetNormalID {
    input:
      bam = normal_bam,
      bam_index = normal_bam_index,
      max_retries = standard_runtime.max_retries,
      preemptible = standard_runtime.preemptible,
      mem = standard_runtime.mem,
      cpu = standard_runtime.cpu,
      boot_disk_size = standard_runtime.boot_disk_size,
      disk_pad = standard_runtime.initial_disk_size
  }

  call ContEst {
    input:
      pair_name = pair_name,
      tumor_bam = tumor_bam,
      tumor_bam_index = tumor_bam_index,
      normal_bam = normal_bam,
      normal_bam_index = normal_bam_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      target_intervals = target_intervals,
      snp6_bed = snp6_bed,
      hapmap_vcf = hapmap_vcf,
      disk_size = ref_size + tumor_size + normal_size,
      runtime_params = standard_runtime
  }

  call deTiN {
    input:
      pair_name = pair_name,
      tumor_name = GetTumorID.sample_id,
      normal_name = GetNormalID.sample_id,
      callstats_file = input_callstats,
      vcf_file = input_vcf,
      seg_file = input_seg,
      tumor_hets = tumor_hets,
      normal_hets = normal_hets,
      exac_pickle = exac_pickle,
      mutation_prior = mutation_prior,
      TiN_prior = TiN_prior,
      runtime_params = standard_runtime
  }

  call OrientationBiasQ as OxoG_Q {
    input:
      stub = "oxog",
      pair_name = pair_name,
      pre_adapter_detail_metrics = pre_adapter_detail_metrics,
      runtime_params = standard_runtime
    }

  call OrientationBiasQ as FFPE_Q {
    input:
      stub = "ffpe",
      pair_name = pair_name,
      pre_adapter_detail_metrics = pre_adapter_detail_metrics,
      runtime_params = standard_runtime
  }


  output {
    # Cross-Sample Contamination
    File contamination_data=ContEst.contest_output
    File contest_base_report = ContEst.contest_base_report
    Float frac_contam = ContEst.contest_frac_contam
    String frac_contam_CI = ContEst.contest_frac_contam_CI

    # deTiN outputs
    Float TiN = deTiN.TiN
    Int deTiN_rescued_SSNVs = deTiN.deTiN_rescued_SSNVs
    String TiN_CI = deTiN.TiN_CI
    File deTiN_call_stats = deTiN.deTiN_call_stats
    File deTiN_SSNVs_plot = deTiN.deTiN_SSNVs_plot
    File deTiN_aSCNA_model = deTiN.deTiN_aSCNA_model
    File deTiN_aSCNA_kmeans_RSS_plot = deTiN.deTiN_aSCNA_kmeans_RSS_plot
    File deTiN_aSCNA_scatter_plot = deTiN.deTiN_aSCNA_scatter_plot
    File deTiN_TiN_modes_plot = deTiN.deTiN_TiN_modes_plot
    File deTiN_indels = deTiN.deTiN_indels
    File deTiN_segments = deTiN.deTiN_segments

    # Orientation bias Q-values
    Float oxog_q_value = OxoG_Q.q_value
    Float ffpe_q_value = FFPE_Q.q_value
  }
}

task ContEst {
  input {
    String pair_name
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File target_intervals    
    File snp6_bed
    File hapmap_vcf

    Int disk_size
    Runtime runtime_params
  }

  Int compute_mem = runtime_params.mem * 1000 - 500
  Int disk_gb = runtime_params.initial_disk_size + disk_size +
                ceil(size(target_intervals,"GB") + size(snp6_bed,"GB") + size(hapmap_vcf,"GB"))

  command <<<
    set -euxo pipefail

    java "-Xmx~{compute_mem}m" -Djava.io.tmpdir=/tmp -jar /root/GATK_v3.5-0-g36282e4.jar \
      -T ContEst \
      -I:eval ~{tumor_bam} \
      -I:genotype ~{normal_bam} \
      -L ~{target_intervals} \
      -L ~{snp6_bed} \
      -isr INTERSECTION \
      -R ~{ref_fasta} \
      -l INFO \
      -pf ~{hapmap_vcf} \
      --trim_fraction 0.03 \
      --beta_threshold 0.05 \
      -o ~{pair_name}.ContEst.out.txt \
      -br ~{pair_name}.ContEst.base_report.txt \
      -mbc 100 \
      --min_genotype_depth 30 \
      --min_genotype_ratio 0.8

    python <<CODE
    import csv

    def frac(val):
        try:
            return str(float(val)/100.0)
        except ValueError:
            return "-1"

    fi = open('~{pair_name}.ContEst.out.txt')
    d = [elem for elem in csv.DictReader(fi, delimiter='\t')][0]

    # convert percentage to fraction, and output the contam frac and the CI
    with open('frac_contam.txt', 'w') as fa:
        fa.write( frac(d['contamination']) )
    with open('frac_contam_CI.txt', 'w') as fb:
        fb.write( frac(d['confidence_interval_95_low']) + ' - ' +
                  frac(d['confidence_interval_95_high']) )
    CODE
  >>>

  output {
    File contest_output = "~{pair_name}.ContEst.out.txt"
    File contest_base_report = "~{pair_name}.ContEst.base_report.txt"
    Float contest_frac_contam = read_float("frac_contam.txt")
    String contest_frac_contam_CI = read_string("frac_contam_CI.txt")
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


task deTiN {
  input {
    String pair_name
    String tumor_name
    String normal_name

    File? callstats_file
    File vcf_file
    File seg_file
    File tumor_hets
    File normal_hets

    File exac_pickle
    Float mutation_prior
    Float TiN_prior

    Runtime runtime_params
  }

  Int disk_gb = runtime_params.initial_disk_size + ceil(size(tumor_hets,"GB") + size(normal_hets,"GB")
                + size(exac_pickle,"GB") + size(seg_file,"GB") + size(vcf_file,"GB") + size(callstats_file,"GB"))
  Boolean is_vcf_compressed = basename(vcf_file, ".gz") != basename(vcf_file)

  command <<<
    set -euxo pipefail

    # generate null outputs
    echo null > ~{pair_name}.TiN_estimate.txt
    echo null > ~{pair_name}.TiN_estimate_CI.txt
    echo null > ~{pair_name}.deTiN_SSNVs.txt
    echo null > ~{pair_name}.deTiN_indels.txt
    echo null > ~{pair_name}.deTiN_aSCNAs.txt
    echo null > ~{pair_name}_SSNVs_plot.png
    echo null > ~{pair_name}_KmeansEval_plot.png
    echo null > ~{pair_name}_TiN_models_plot.png
    echo null > ~{pair_name}_TiN_hets_aSCNA_model.png
    echo null > ~{pair_name}_KmeansEval_scatter_plot.png

    # decompress VCF and create deTiN input files
    INPUT_VCF="~{pair_name}.vcf"
    SNV_FILE="~{tumor_name}.SNV.callstats.txt"
    INDEL_FILE="~{tumor_name}.INDEL.vcf"
    INDEL_TYPE="~{tumor_name}.INDEL_TYPE.txt"
    if [[ ~{is_vcf_compressed} ]] ; then
      gunzip -c ~{vcf_file} > $INPUT_VCF
    else
      cp ~{vcf_file} $INPUT_VCF
    fi
    python /root/scripts/detin_inputs.py ~{tumor_name} ~{normal_name} $INPUT_VCF
    INDEL_TYPE=`cat $INDEL_TYPE`

    # overwrite SNV input file if there is an input callstats file
    if [[ -f "~{callstats_file}" ]] ; then
      cp ~{callstats_file} $SNV_FILE      
    fi

    python /root/deTiN-2.0.1/deTiN/deTiN.py \
        --mutation_data_path $SNV_FILE \
        --cn_data_path ~{seg_file} \
        --tumor_het_data ~{tumor_hets} \
        --normal_het_data ~{normal_hets} \
        --exac_data_path ~{exac_pickle} \
        --output_name ~{pair_name} \
        --TiN_prior ~{TiN_prior} \
        --mutation_prior ~{mutation_prior} \
        --indel_data_path $INDEL_FILE \
        --indel_data_type $INDEL_TYPE
  >>>

  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }

  output {
    Float TiN = read_float("~{pair_name}.TiN_estimate.txt")
    String TiN_CI = read_string("~{pair_name}.TiN_estimate_CI.txt")
    Int deTiN_rescued_SSNVs = read_int("~{pair_name}.number_of_SSNVs_added.txt")
    File deTiN_call_stats = "~{pair_name}.deTiN_SSNVs.txt"
    File deTiN_SSNVs_plot = "~{pair_name}_SSNVs_plot.png"
    File deTiN_aSCNA_kmeans_RSS_plot = "~{pair_name}_KmeansEval_plot.png"
    File deTiN_aSCNA_scatter_plot = "~{pair_name}_KmeansEval_scatter_plot.png"
    File deTiN_TiN_modes_plot = "~{pair_name}_TiN_models_plot.png"
    File deTiN_indels = "~{pair_name}.deTiN_indels.txt"
    File deTiN_aSCNA_model = "~{pair_name}_TiN_hets_aSCNA_model.png"
    File deTiN_segments = "~{pair_name}.deTiN_aSCNAs.txt"
  }

  meta {
    author: "Junko Tsuji"
  }
}

task OrientationBiasQ {
  input {
    String pair_name
    String stub
    File pre_adapter_detail_metrics
    Runtime runtime_params
  }

  Int disk_gb = runtime_params.initial_disk_size + ceil(size(pre_adapter_detail_metrics,"GB"))

  command <<<
    set -e

    # determine orientation bias contexts
    if [[ "~{stub}" = "oxog" ]] ; then
      echo "Artifact mode: OxoG"
      CONTEXT=".GG"
      REF_ALLELE="G"
      REF_ALLELE_COMP="C"
      ARTIFACT_ALLELE="T"
      ARTIFACT_ALLELE_COMP="A"
    elif [[ "~{stub}" = "ffpe" ]] ; then
      echo "Artifact mode: FFPE"
      CONTEXT=".CG"
      REF_ALLELE="C"
      REF_ALLELE_COMP="G"
      ARTIFACT_ALLELE="T"
      ARTIFACT_ALLELE_COMP="A"
    else
      echo "Error: unknown artifact mode" 1>&2
      exit 1
    fi

    # parse pre_adapter_detail_metrics for extracting Q-value
    python /root/scripts/annotate_orientationBiasQ.py \
      -i ~{pair_name} \
      -m ~{pre_adapter_detail_metrics} \
      -c $CONTEXT \
      -a $ARTIFACT_ALLELE
  >>>

  output {
    Float q_value=read_float("${pair_name}.orientation_BiasQ.txt")
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
