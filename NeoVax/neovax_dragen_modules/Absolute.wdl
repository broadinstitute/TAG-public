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

workflow Absolute {
  input {
    String pair_name
    File seg_file
    File maf_file
    Float skew

    # workflow level runtime parameters
    Int preemptible = 2
    Int max_retries = 1
    Int additional_disk = 20
    Int boot_disk_size = 15
    Int mem = 8
    Int cpu = 1
  }

  Runtime standard_runtime = { "preemptible": preemptible,
                                "max_retries": max_retries,
                                "mem": mem,
                                "cpu": cpu,
                                "docker": "us.gcr.io/tag-public/neovax-tag-cga:v1",
                                "boot_disk_size": boot_disk_size,
                                "initial_disk_size": additional_disk }

  call RunAbsolute {
    input:
      analysis_id = pair_name,
      seg_file = seg_file,
      maf_file = maf_file,
      skew = skew,
      runtime_params = standard_runtime
  }
  
  output {
    File absolute_highres_plot = RunAbsolute.absolute_highres_plot
    File absolute_rdata = RunAbsolute.absolute_rdata
  }
}

task RunAbsolute {
  input {
    String analysis_id
    File seg_file
    File maf_file
    Float skew
    Runtime runtime_params
  }

  Int disk_gb = runtime_params.initial_disk_size + ceil(size(seg_file, "GB") + size(maf_file, "G"))

  command <<<
    set -euxo pipefail

    # correct [Start|End]_Position to [Start|End]_position
    grep ^Hugo_Symbol ~{maf_file} | tr "\t" "\n" > ~{analysis_id}.column_names.txt
    python /usr/local/bin/subset_maf_columns.py --correct-maf-dialect \
      ~{analysis_id}.column_names.txt ~{maf_file} ~{analysis_id}.corrected.maf

    SNV_MAF="~{analysis_id}.snv.maf"
    INDEL_MAF="~{analysis_id}.indel.maf"
    python /usr/local/bin/split_maf_indel_snp.py -i ~{analysis_id}.corrected.maf -o $SNV_MAF -f Variant_Type -v "SNP|DNP|TNP|MNP|ONP"
    python /usr/local/bin/split_maf_indel_snp.py -i ~{analysis_id}.corrected.maf -o $INDEL_MAF -f Variant_Type -v "INS|DEL"
        
    grep -v "NA" ~{seg_file} > no_nan_segs.tsv

    Rscript /xchip/tcga/Tools/absolute/releases/v1.5/run/ABSOLUTE_cli_start.R \
      --seg_dat_fn no_nan_segs.tsv \
      --maf_fn $SNV_MAF \
      --indelmaf_fn $INDEL_MAF \
      --sample_name ~{analysis_id} \
      --results_dir . \
      --ssnv_skew ~{skew} \
      --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5/
  >>>
  
  output {
    File absolute_highres_plot = "${analysis_id}.ABSOLUTE_plot.pdf"
    File absolute_rdata = "${analysis_id}.PP-modes.data.RData"
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
    email: "jtsuji@broadinstitute.org"
    description: "Adopted from [the CGA WES Characterization Workflow](https://firecloud.terra.bio/#workspaces/broad-fc-getzlab-workflows/CGA_WES_Characterization_OpenAccess) and tweaked for handling ONPs and MNPs"
  }
}
