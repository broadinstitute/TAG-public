version 1.0

workflow QUICviz {
    input {
        String sampleID
        String tumorType
        String quicvizDocker = "us-central1-docker.pkg.dev/tag-team-160914/gptag-dockers/cmi_quicviz:0.4.1"
        File allelicCountsNormal
        File allelicCountsTumor
        File denoisedCopyRatiosNormal
        File denoisedCopyRatiosTumor
        File calledCopyRatioSegTumor
        File oncotatedCalledTumor
    }
    call QUICviz {
        input:
            sampleID = sampleID,
            tumorType = tumorType,
            quicvizDocker = quicvizDocker,
            allelicCountsNormal = allelicCountsNormal,
            allelicCountsTumor = allelicCountsTumor,
            denoisedCopyRatiosNormal = denoisedCopyRatiosNormal,
            denoisedCopyRatiosTumor = denoisedCopyRatiosTumor,
            calledCopyRatioSegTumor = calledCopyRatioSegTumor,
            oncotatedCalledTumor = oncotatedCalledTumor
    }
    output {
        File QUICvizPDF = QUICviz.QUICvizPDF
        File GeneLevelCNV = QUICviz.GeneLevelCNV
        File AllChrPlot = QUICviz.AllChrPlot
    }

    meta {
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        description: "QUICviz.wdl is based on the QUICviz_v0.4 R script developed by Alex Neil, which is a tool for visualizing CNV data"
    }
}

task QUICviz {
    input {
        String sampleID
        String tumorType
        String quicvizDocker
        File allelicCountsNormal
        File allelicCountsTumor
        File denoisedCopyRatiosNormal
        File denoisedCopyRatiosTumor
        File calledCopyRatioSegTumor
        File oncotatedCalledTumor
        Int memory = 16
        Int cpu = 4
    }
    command <<<
        set -e
        mkdir outputs

        Rscript /BaseImage/CMI_QUICviz/scripts/CMI_QUICviz_v0.3.R \
            --sample ~{sampleID} \
            --tumor_type ~{tumorType} \
            --normal_acf ~{allelicCountsNormal} \
            --normal_cr ~{denoisedCopyRatiosNormal} \
            --tumor_acf ~{allelicCountsTumor} \
            --tumor_cr ~{denoisedCopyRatiosTumor} \
            --tumor_cr_seg ~{calledCopyRatioSegTumor} \
            --tumor_seg_oncotated ~{oncotatedCalledTumor} \
            --output_dir outputs/


    >>>
    output {
        File QUICvizPDF = "outputs/chromosome_plots.pdf"
        File GeneLevelCNV = "outputs/gene_level_calls.csv"
        File AllChrPlot = "outputs/All_chr.png"
    }
    runtime {
        docker: quicvizDocker
        memory: memory + " GB"
        cpu: cpu
        disks: "local-disk 100 HDD"
    }
}
