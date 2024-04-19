version 1.0

workflow QUICviz {
    input {
        String sampleID
        String tumorType
        String quicvizDocker = "us-central1-docker.pkg.dev/tag-team-160914/gptag-dockers/cmi_quicviz:0.3.0"
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
    call mergeImages {
        input:
            SampleID = sampleID,
            TumorType = tumorType,
            plot = QUICviz.plot
    }
    output {
        File QUICvizPDF = mergeImages.compiled_pdf
    }
    meta {
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        description: "QUICviz.wdl is based on the QUICviz_v0.3 R script developed by Alex Neil, which is a tool for visualizing CNV data"
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
    }
    command <<<
        set -e
        mkdir -p /BaseImage/CMI_QUICviz/outputs/

        Rscript /BaseImage/CMI_QUICviz/scripts/CMI_QUICviz_v0.3.R \
            --sample ~{sampleID} \
            --tumor_type ~{tumorType} \
            --normal_acf ~{allelicCountsNormal} \
            --normal_cr ~{denoisedCopyRatiosNormal} \
            --tumor_acf ~{allelicCountsTumor} \
            --tumor_cr ~{denoisedCopyRatiosTumor} \
            --tumor_cr_seg ~{calledCopyRatioSegTumor} \
            --tumor_seg_oncotated ~{oncotatedCalledTumor} \
            --output_dir /BaseImage/CMI_QUICviz/outputs/

    >>>
    output {
        Array[File] plot = glob("outputs/*.png")
    }
    runtime {
        docker: quicvizDocker
        memory: "16 GB"
        cpu: "4"
        disks: "local-disk 100 HDD"
    }
}
task mergeImages {
    input {
        String SampleID
        String TumorType
        Array[File] plot
    }
    command <<<
        for i in ~{plot}; do mv $i /output/images/; done

        python <<CODE
        pip3 install img2pdf
        import img2pdf
        import glob
        import os

        with open(f"/output/~{SampleID}_~{TumorType}_QUICviz.pdf","wb") as f:
            f.write(img2pdf.convert(glob.glob("/output/images/*.jpg")))

        CODE
    >>>
    output {
        File compiled_pdf = "/output/~{SampleID}_~{TumorType}_QUICviz.pdf"
    }
    runtime {
        docker: "us.gcr.io/tag-team-160914/neovax-parsley:2.2.1.0"
        memory: "16 GB"
        cpu: "4"
        disks: "local-disk 100 HDD"
    }
}