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
            plot = QUICvizPlots,
            quicvizDocker = quicvizDocker
    }
    output {
        File QUICvizPDF = mergeImages.chr_pdf
        File AllChrPlot = mergeImages.allchr_plot
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
        Array[File] plot = glob("outputs/*.png")
    }
    runtime {
        docker: quicvizDocker
        memory: memory
        cpu: cpu
        disks: "local-disk 100 HDD"
    }
}
task mergeImages {
    input {
        String SampleID
        String TumorType
        Array[File] plot
        String quicvizDocker
        Int memory = 16
        Int cpu = 4
    }
    command <<<
        mkdir -p output/images
        for i in `ls ~{sep=" " plot}`; do mv $i output/images/; done

        python3 <<CODE
        import img2pdf
        import glob
        import os

        # Get list of PNG files sorted
        png_files = sorted(glob.glob("output/images/*.png"))
        numeric_png_files = [file for file in png_files if os.path.basename(file).split('.')[0].isdigit()]
        png_files = sorted(numeric_png_files, key=lambda x: int(os.path.basename(x).split('.')[0]))

        with open(f"output/~{SampleID}_~{TumorType}_QUICviz.pdf","wb") as f:
            f.write(img2pdf.convert(png_files))
        CODE
    >>>
    output {
        File chr_pdf = "output/~{SampleID}_~{TumorType}_QUICviz.pdf"
        File allchr_plot = "output/images/All_chr.png"
    }
    runtime {
        docker: quicvizDocker
        memory: memory
        cpu: cpu
        disks: "local-disk 100 HDD"
    }
}