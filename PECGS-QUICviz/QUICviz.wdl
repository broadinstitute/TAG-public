version 1.0

workflow QUICviz {
    input {
        String sampleID
        Boolean isPECGS = true
        String tumorType
        String quicvizDocker = "us-central1-docker.pkg.dev/tag-team-160914/gptag-dockers/cmi_quicviz:0.4.2"
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
        Boolean isPECGS
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
        Int maxRetries = 3
    }
    command <<<
        set -e
        mkdir outputs

        if ~{isPECGS}; then
            IFS='-' read -r tumor_sample normal_sample <<< "~{sampleID}"
            echo "Input Tumor Sample: $tumor_sample"
            echo "Input Normal Sample: $normal_sample"
        else
            tumor_sample=~{sampleID}
            echo "Input Tumor Sample: $tumor_sample"
        fi


        Rscript /BaseImage/CMI_QUICviz/scripts/CMI_QUICviz_v0.4.1.R \
            --sample $tumor_sample \
            --tumor_type ~{tumorType} \
            --normal_acf ~{allelicCountsNormal} \
            --normal_cr ~{denoisedCopyRatiosNormal} \
            --tumor_acf ~{allelicCountsTumor} \
            --tumor_cr ~{denoisedCopyRatiosTumor} \
            --tumor_cr_seg ~{calledCopyRatioSegTumor} \
            --tumor_seg_oncotated ~{oncotatedCalledTumor} \
            --output_dir outputs/

        mv outputs/*chromosome_plots.pdf outputs/~{sampleID}_chromosome_plots.pdf
        mv outputs/*gene_level_calls.csv outputs/~{sampleID}_gene_level_calls.csv
        mv outputs/*_all_chr.png outputs/~{sampleID}_All_chr.png

    >>>
    output {
        File QUICvizPDF = "outputs/~{sampleID}_chromosome_plots.pdf"
        File GeneLevelCNV = "outputs/~{sampleID}_gene_level_calls.csv"
        File AllChrPlot = "outputs/~{sampleID}_All_chr.png"
    }
    runtime {
        docker: quicvizDocker
        memory: memory + " GB"
        cpu: cpu
        disks: "local-disk 100 HDD"
        maxRetries: maxRetries
    }
}
