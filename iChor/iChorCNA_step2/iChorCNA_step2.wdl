task bundle_plots {
    Array[File] CNA_plots
    File summary_script
    String sample_set_name

    Int? disk_size

    command {
        python ${summary_script} --prefix ${sample_set_name} ${sep=' ' CNA_plots}
        pdflatex ${sample_set_name}.tex
    }
    runtime {
        docker: "us.gcr.io/tag-team-160914/tag-tools:0.0.4"
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
