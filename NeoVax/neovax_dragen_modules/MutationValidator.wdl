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

workflow MutationValidator {
  input {
    String pair_name
    

    # type of sequencing method used for MAF
    # options: 'wgs', 'wex', 'rna', 'targeted', 'low pass'
    String maf_type = "wex"
    File maf_file
    
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index
    File tumor_rna_bam
    File tumor_rna_bam_index
    Boolean add_chr_prefix_to_rna = false

    # workflow level runtime parameters
    Int preemptible = 2
    Int max_retries = 1
    Int additional_disk = 80
    Int boot_disk_size = 15
    Int mem = 10
    Int cpu = 1
  }    

  Runtime standard_runtime = { "preemptible": preemptible,
                                "max_retries": max_retries,
                                "mem": mem,
                                "cpu": cpu,
                                "docker": "us.gcr.io/tag-public/neovax-tag-cga:v1",
                                "boot_disk_size": boot_disk_size,
                                "initial_disk_size": additional_disk }

  call RunMutationValidator {
    input:
      pair_name = pair_name,
      maf_file = maf_file,
      maf_type = maf_type,
      add_chr_prefix_to_rna = add_chr_prefix_to_rna,
      tumor_bam = tumor_bam,
      tumor_bam_index = tumor_bam_index,
      normal_bam = normal_bam,
      normal_bam_index = normal_bam_index,
      tumor_rna_bam = tumor_rna_bam,
      tumor_rna_bam_index = tumor_rna_bam_index,
      runtime_params = standard_runtime
  }
  output{
    File validated_maf = RunMutationValidator.validated_maf
    File pileup_preprocessing_txt = RunMutationValidator.pileup_preprocessing_txt
  }
}

task RunMutationValidator {
  input {
    String pair_name
    File maf_file
    String maf_type
    Boolean add_chr_prefix_to_rna
    File tumor_bam
    File tumor_bam_index
    File normal_bam
    File normal_bam_index
    File tumor_rna_bam
    File tumor_rna_bam_index
    Runtime runtime_params
  }

  Int tumor_size = ceil(size(tumor_bam,"GB") + size(tumor_bam_index,"GB"))
  Int normal_size = ceil(size(normal_bam,"GB") + size(normal_bam_index,"GB"))
  Int rna_size = ceil(size(tumor_rna_bam,"GB") + size(tumor_rna_bam_index,"GB"))

  Int disk_gb = runtime_params.initial_disk_size + tumor_size + normal_size + rna_size +
                ceil(size(maf_file,"GB"))
  
  command <<<
    set -euxo pipefail

    # correct [Start|End]_Position to [Start|End]_position
    grep ^Hugo_Symbol ~{maf_file} | tr "\t" "\n" > ~{pair_name}.column_names.txt
    python /usr/local/bin/subset_maf_columns.py --correct-maf-dialect \
      ~{pair_name}.column_names.txt ~{maf_file} ~{pair_name}.corrected.maf

    SNV_ONLY_MAF="~{pair_name}.snv_only.maf"
    MNP_MAF="~{pair_name}.mnp.maf"
    MNP_SPLIT_MAF="~{pair_name}.mnp_split.maf"
    
    SNV_MAF="~{pair_name}.snv.maf"
    INDEL_MAF="~{pair_name}.indel.maf"
    
    # split MAF into INDELs, SNPs, and MNPs/ONPs
    python /usr/local/bin/split_maf_indel_snp.py -i ~{pair_name}.corrected.maf -o $INDEL_MAF -f Variant_Type -v "INS|DEL"
    python /usr/local/bin/split_maf_indel_snp.py -i ~{pair_name}.corrected.maf -o $SNV_ONLY_MAF -f Variant_Type -v "SNP"
    python /usr/local/bin/split_maf_indel_snp.py -i ~{pair_name}.corrected.maf -o $MNP_MAF -f Variant_Type -v "DNP|TNP|MNP|ONP"

    # separate MNPs/ONPs into SNPs and merge the MNP-based SNP MAF with SNP-only MAF
    { MNP_NUMBER=`grep -v ^"#" $MNP_MAF | grep -v ^Hugo_Symbol | wc -l` || :; }
    if [[ $MNP_NUMBER -gt 0 ]] ; then
      python /usr/local/bin/split_onp_to_snp.py $MNP_MAF $MNP_SPLIT_MAF
      python /usr/local/bin/maf_merge.py $SNV_ONLY_MAF $MNP_SPLIT_MAF $SNV_MAF
    else
      cp $SNV_ONLY_MAF $SNV_MAF
    fi
    
    # generate a pileup file as a processed file
    RNATYPE="hg19"
    if [[ ~{add_chr_prefix_to_rna} = true ]]; then
      RNATYPE="hg19-chr"
    fi
    python /opt/src/fh_MutationValidatorPreprocess/validation_wrapper_firehose_library_hack.py \
      --mafsnp $SNV_MAF \
      --mafindel $INDEL_MAF \
      --wextumor ~{tumor_bam} \
      --wexnormal ~{normal_bam} \
      --rnatumor ~{tumor_rna_bam} \
      --rnatype $RNATYPE \
      --out ~{pair_name}

    PREPROCESSED_FILE="~{pair_name}.pileup_preprocessing.txt"
    VALIDATED_SNV_MAF="~{pair_name}.snv.validated.maf"
    VALIDATED_MNP_MAF="~{pair_name}.mnp.validated.maf"
    VALIDATED_INDEL_MAF="~{pair_name}.indel.validated.maf"
    VALIDATED_MAF="~{pair_name}.validated.maf"

    bash -c "source /matlab_source_file_2012a.sh && /opt/src/fh_MutationValidator/mutation_validator_wrapper $PREPROCESSED_FILE $SNV_MAF ~{maf_type} ~{pair_name}.snv"
    bash -c "source /matlab_source_file_2012a.sh && /opt/src/fh_MutationValidator/mutation_validator_wrapper $PREPROCESSED_FILE $INDEL_MAF ~{maf_type} ~{pair_name}.indel"

    # integrate the validated SNV information into MNPs
    if [[ $MNP_NUMBER -gt 0 ]] ; then
      python /usr/local/bin/chain_snp_to_onp.py --columns-to-chain /usr/local/bin/maf_columns/mutation_validator_maf_columns.txt $MNP_MAF $VALIDATED_SNV_MAF $VALIDATED_MNP_MAF
    else
      cp $VALIDATED_SNV_MAF $VALIDATED_MNP_MAF
    fi
    python /usr/local/bin/tsvConcatFiles.py $VALIDATED_MNP_MAF $VALIDATED_INDEL_MAF --outputFilename=$VALIDATED_MAF
  >>>

  output {    
    File pileup_preprocessing_txt="~{pair_name}.pileup_preprocessing.txt"
    File validated_maf="~{pair_name}.validated.maf"
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
    description: "Adopted from [the CGA WES Characterization Workflow](https://firecloud.terra.bio/#workspaces/broad-fc-getzlab-workflows/CGA_WES_Characterization_OpenAccess) and tweaked for handling ONPs and MNPs and different header names"
  }
}
