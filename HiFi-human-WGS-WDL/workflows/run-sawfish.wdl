version 1.0

import "./wdl-common/wdl/tasks/sawfish.wdl" as Sawfish
import "./wdl-common/wdl/structs.wdl"
import "./wdl-common/wdl/workflows/backend_configuration/backend_configuration.wdl" as BackendConfiguration

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
  }

  call BackendConfiguration.backend_configuration {
    input:
      backend            = "GCP",
      container_registry = "quay.io/pacbio"
  }

  RuntimeAttributes default_runtime_attributes = backend_configuration.on_demand_runtime_attributes

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
    # per sample sv signatures
    File discover_tar = sawfish_discover.discover_tar

    # sawfish outputs for single sample
    File? sv_vcf       = sawfish_call.vcf
    File? sv_vcf_index = sawfish_call.vcf_index
    File? supporting_reads = sawfish_call.supporting_reads
  }
}
