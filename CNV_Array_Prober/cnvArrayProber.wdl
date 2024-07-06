version 1.0

workflow cnvArrayProber {
    input{
        String sampleName
        File cnvBedFile
        File CytoSNP850K_Support_Csv
        File GDA_Support_Csv
        String cnvProberDocker = "us.gcr.io/tag-public/cnv-array-prober:0.0.0"
    }
    call cnvArrayProber {
        input:
            sampleName = sampleName,
            cnvBedFile = cnvBedFile,
            CytoSNP850K_Support_Csv = CytoSNP850K_Support_Csv,
            GDA_Support_Csv = GDA_Support_Csv,
            cnvProberDocker = cnvProberDocker
    }
    output{
        File cnvProbeAnnotation = cnvArrayProber.cnvProbeAnnotation
        File cnvProbePlots = cnvArrayProber.cnvProbePlots
    }
    meta {
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        description: "This workflow takes a CNV bed file and CytoSNP-850K and GDA support files as input and outputs a csv file with probe information for each CNV interval. Additionally,  output a PDF file with plots for each CNV interval the number of probes in the CytoSNP-850K and GDA arrays."
    }
}

task cnvArrayProber {
    input{
        String sampleName
        File cnvBedFile
        File CytoSNP850K_Support_Csv
        File GDA_Support_Csv
        String cnvProberDocker
        Int memory = 32
        Int cpu = 2
        Int disk_size_gb = 500
        Boolean use_ssd = false
    }
    command <<<
    set -e
    mkdir output

    conda run --no-capture-output \
            -n prober_env \
            python3 /BaseImage/cnvArrayProber/scripts/cnvArrayProber.py \
            -b ~{cnvBedFile} \
            -c ~{CytoSNP850K_Support_Csv} \
            -g ~{GDA_Support_Csv} \
            -o output/~{sampleName}
    >>>
    output{
        File cnvProbeAnnotation = "output/~{sampleName}CNV_Probe_Mappings.csv"
        File cnvProbePlots = "output/~{sampleName}CNV_Probe_Mappings_Plots.pdf"
    }
    runtime {
        docker: cnvProberDocker
        memory: memory
        cpu: cpu
        disks: "local-disk " + disk_size_gb + if use_ssd then " SSD" else " HDD"
    }
}