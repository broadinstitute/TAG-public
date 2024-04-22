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

    Array[File] QUICvizPlots = QUICviz.plot
    call mergeImages {
        input:
            SampleID = sampleID,
            TumorType = tumorType,
            plot = QUICvizPlots
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

        ls outputs/
        readlink -f outputs/*png
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
        mkdir -p output/images
        for i in `ls ~{sep=" " plot}`; do mv $i output/images/; done
        ls output/images/
        source activate NeoVax-Input-Parser
        pip3 install img2pdf

        python <<CODE
        import img2pdf
        import glob
        import os

        # Get list of PNG files sorted
        png_files = sorted(glob.glob("output/images/*.png"), key=lambda x: int(os.path.basename(x).split('.')[0]))

        with open(f"output/~{SampleID}_~{TumorType}_QUICviz.pdf","wb") as f:
            f.write(img2pdf.convert(png_files))
        CODE
    >>>
    output {
        File compiled_pdf = "output/~{SampleID}_~{TumorType}_QUICviz.pdf"
    }
    runtime {
        docker: "us.gcr.io/tag-team-160914/neovax-parsley:2.2.1.0"
        memory: "16 GB"
        cpu: "4"
        disks: "local-disk 100 HDD"
    }
}