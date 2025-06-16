version 1.0

import "./wdl-common/wdl/tasks/sawfish.wdl" as Sawfish
import "./wdl-common/wdl/structs.wdl"

workflow run_sawfish {
  meta {
    description: "Run sawfish on individual samples."
  }

  input {
    String sample_id
    String sex
    File aligned_bam_data
    File aligned_bam_index

    File ref_map_file

    Boolean single_sample = true

    Boolean gpu

    RuntimeAttributes default_runtime_attributes
  }

  Map[String, String] ref_map = read_map(ref_map_file)


  call Sawfish.sawfish_discover {
    input:
      sex                 = sex,
      aligned_bam         = aligned_bam_data,
      aligned_bam_index   = aligned_bam_index,
      ref_fasta           = ref_map["fasta"],                           # !FileCoercion
      ref_index           = ref_map["fasta_index"],                     # !FileCoercion
      out_prefix          = "~{sample_id}.~{ref_map['name']}",
      expected_male_bed   = ref_map["hificnv_expected_bed_male"],       # !FileCoercion
      expected_female_bed = ref_map["hificnv_expected_bed_female"],     # !FileCoercion
      runtime_attributes  = default_runtime_attributes
  }

  if (single_sample) {
    call Sawfish.sawfish_call {
      input:
        discover_tars       = [sawfish_discover.discover_tar],
        aligned_bams        = [aligned_bam_data],
        aligned_bam_indices = [aligned_bam_index],
        ref_fasta           = ref_map["fasta"],                                      # !FileCoercion
        ref_index           = ref_map["fasta_index"],                                # !FileCoercion
        out_prefix          = "~{sample_id}.~{ref_map['name']}.structural_variants",
        runtime_attributes  = default_runtime_attributes
    }
  }

  output {
    # mosdepth outputs
    File   mosdepth_summary                 = mosdepth.summary
    File   mosdepth_region_bed              = mosdepth.region_bed
    File   mosdepth_region_bed_index        = mosdepth.region_bed_index
    File   mosdepth_depth_distribution_plot = mosdepth.depth_distribution_plot
    String inferred_sex                     = mosdepth.inferred_sex
    String stat_mean_depth                  = mosdepth.stat_mean_depth

    # per sample sv signatures
    File discover_tar = sawfish_discover.discover_tar

    # sawfish outputs for single sample
    File? sv_vcf       = sawfish_call.vcf
    File? sv_vcf_index = sawfish_call.vcf_index
    File? supporting_reads = sawfish_call.supporting_reads
  }
}
