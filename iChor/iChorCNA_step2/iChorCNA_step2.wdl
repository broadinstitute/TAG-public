task bundle_plots {
    Array[File] CNA_plots
    File summary_script
    String sample_set_name
    String docker = "us.gcr.io/tag-public/tag-tools:1.0.0"
    Int? disk_size

    command {
        python ${summary_script} --prefix ${sample_set_name} ${sep=' ' CNA_plots}
        pdflatex ${sample_set_name}.tex
    }
    runtime {
        docker: docker
        memory: "2 GB"
        cpu: 1
        disks: "local-disk " + select_first([disk_size, 20]) + " HDD"
    }
    output {
        File bundledPDF = "${sample_set_name}.pdf"
    }
}

workflow ichorSummaryBundle {
    Array[File] CNA_plots
    File summary_script
    String sample_set_name

    call bundle_plots {
        input: CNA_plots = CNA_plots,
               summary_script = summary_script,
               sample_set_name = sample_set_name
    }
    output {
        bundle_plots.*
    }
}
