version 1.0

workflow runTPMCutoff {
  call findTPMCutoff

  output{
    File cutoff_plot_pdf = findTPMCutoff.cutoff_plot_pdf
    File cutoff_plot_png = findTPMCutoff.cutoff_plot_png
    Float first_cutoff_tpm = findTPMCutoff.first_cutoff_tpm
    Int first_cutoff_genes = findTPMCutoff.first_cutoff_genes
    Float second_cutoff_tpm = findTPMCutoff.second_cutoff_tpm
    Int second_cutoff_genes = findTPMCutoff.second_cutoff_genes
  }
}

task findTPMCutoff {
  input{
  File tpm_cutoff_script
  String patient_id
  String disease_code
  File genes_tpm

  String? docker_override
  Int? extra_disk
  Float? extra_mem
  Int? maxRetries
  Int? preemptible
  }
  String output_basename = "NeoVax_" + disease_code + "_" + patient_id
  String docker = select_first([docker_override, "us.gcr.io/tag-public/tag-tools:1.0.0"])
  
  command {
    Rscript ~{tpm_cutoff_script} ~{genes_tpm} ~{patient_id} ~{disease_code}
  }

  runtime {
    docker: docker
    disk: 50 + select_first([0, extra_disk]) + " GB"
    memory: 5 + select_first([0, extra_mem]) + "GB"
    maxRetries: select_first([2, maxRetries])
    preemptible: select_first([2, preemptible])
  }

  output {
    File cutoff_plot_pdf = "~{output_basename}.cutoff_plot.pdf"
    File cutoff_plot_png = "~{output_basename}.cutoff_plot.png"
    Float first_cutoff_tpm = read_float("~{output_basename}.first_tpm_cutoff.txt")
    Int first_cutoff_genes = read_int("~{output_basename}.first_tpm_cutoff_genes_detected.txt")
    Float second_cutoff_tpm = read_float("~{output_basename}.second_tpm_cutoff.txt")
    Int second_cutoff_genes = read_int("~{output_basename}.second_tpm_cutoff_genes_detected.txt")
  }
}