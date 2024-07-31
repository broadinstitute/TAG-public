version 1.0

workflow CNV_Profiler {
    input{
        String sampleName
        String cnvProfiler_Docker = "us-central1-docker.pkg.dev/tag-team-160914/gptag-dockers/covprofileviz:0.0.3"
        File cramOrBamFile
        File cramOrBamIndexFile
        File referenceFasta
        File referenceFastaIndex
        File referenceDict
        File cnvBedFile
        Boolean heterozygosityCheck = false
        File? hardFilteredVcfFile
    }
    if (basename(cramOrBamFile) != basename(cramOrBamFile, ".cram")) {
        call CramToBam {
            input:
                sampleName = sampleName,
                cramFile = cramOrBamFile,
                cramIndexFile = cramOrBamIndexFile,
                referenceFasta = referenceFasta,
                referenceFastaIndex = referenceFastaIndex,
                referenceDict = referenceDict
        }
    }
    File alignedBam = select_first([cramOrBamFile, CramToBam.output_bam])
    File alignedBai = select_first([cramOrBamIndexFile, CramToBam.output_bai])
    call GetPaddedCnvBed {
        input:
            cnvBedFile = cnvBedFile,
            referenceDict = referenceDict,
            cnvProfiler_Docker = cnvProfiler_Docker
    }
    call SamtoolsDepth {
        input:
            sampleName = sampleName,
            alignedBam = alignedBam,
            alignedBai = alignedBai,
            target_bed = GetPaddedCnvBed.paddedCnvBed

    }
    call cnvDepthProfiler {
        input:
            sampleName = sampleName,
            depthProfile = SamtoolsDepth.depth_profile,
            cnvBedFile = cnvBedFile,
            cnvProfiler_Docker = cnvProfiler_Docker,
            PaddedcnvBedFile = GetPaddedCnvBed.paddedCnvBed
    }
    if (heterozygosityCheck) {
        call HeterozygosityCheck {
            input:
                sampleName = sampleName,
                hardFilteredVcfFile = hardFilteredVcfFile,
                cnvBedFile = cnvBedFile,
                cnvProfiler_Docker = cnvProfiler_Docker
        }
    }
    output {
    File samtools_depth_profile = SamtoolsDepth.depth_profile
    Array[File] cnv_depth_profile = cnvDepthProfiler.cnv_depth_profile
    Array[File]? heterozygosity_plot = HeterozygosityCheck.heterozygosity_plot
    }
    meta {
        description: "This workflow takes a BAM or CRAM file and a CNV bed file as input and generates a coverage profile for the CNV regions in the bed file. Optionally, it can also generate a heterozygosity plot using a hard-filtered VCF file."
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        }
}


task CramToBam {
    input {
      File referenceFasta
      File referenceFastaIndex
      File referenceDict
      #cram and crai must be optional since Normal cram is optional
      File? cramFile
      File? cramIndexFile
      String sampleName
      Int disk_size = 500
      Int mem = 64
      Int cpu = 8
    }

    Int machine_mem = if defined(mem) then mem * 1000 else 6000

    #Calls samtools view to do the conversion
    command <<<
        set -e
        set -o pipefail

        samtools view -h -T ~{referenceFasta} ~{cramFile} |
            samtools view -b -o ~{sampleName}.bam -
        samtools index -b ~{sampleName}.bam
        mv ~{sampleName}.bam.bai ~{sampleName}.bai
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
        cpu: cpu
        memory: machine_mem + " MB"
        disks: "local-disk " + disk_size + " SSD"
    }

    output {
        File output_bam = "~{sampleName}.bam"
        File output_bai = "~{sampleName}.bai"
    }
}

task GetPaddedCnvBed {
    input {
        File cnvBedFile
        File referenceDict
        String cnvProfiler_Docker
        Int mem_gb = 1
        Int cpu = 1
        Int disk_size_gb = 10
    }

    command <<<
        source activate env_viz
        python3 <<CODE
        import re

        # Read the reference dictionary file to get the chromosome lengths
        length_dict = {}
        with open('~{referenceDict}', 'r') as f:
            for i in f.readlines():
                if i.startswith('@SQ\tSN'):
                    chrom = i.split('\t')[1].split('SN:')[1]
                    match = re.match(pattern=r'(chr)?(\d+|X|Y)\b', string=chrom)
                    if match:
                        length = i.split('\t')[2].split('LN:')[1]
                        length_dict[chrom] = int(length)

        # Get CNV interval with padding
        # Padding is 2 times the length of the CNV unless it goes beyond the chromosome length
        padded_cnv_interval_list = []
        with open("~{cnvBedFile}", 'r') as f:
            for line in f:
                chr = line.strip().split('\t')[0]
                start = line.strip().split('\t')[1]
                end = line.strip().split('\t')[2]
                svlen = int(end) - int(start) + 1
                initial_padding_start = int(start)-svlen*2
                initial_padding_end = int(end)+svlen*2
                if initial_padding_start < 0:
                    padding_start = 0
                else:
                    padding_start = initial_padding_start
                if initial_padding_end > length_dict[chr]:
                    padding_end = length_dict[chr]
                else:
                    padding_end = initial_padding_end
                # Add the padded interval to the list
                padded_cnv_interval_list.append(f'{chr}:{padding_start}-{padding_end}')


        with open('padded_cnv.bed', 'a') as f:
            for interval in padded_cnv_interval_list:
                chr = interval.split(':')[0]
                start = interval.split(':')[1].split('-')[0]
                end = interval.split(':')[1].split('-')[1]
                f.write(f"{chr}\t{start}\t{end}" + '\n')
        CODE
    >>>
    runtime {
        docker: cnvProfiler_Docker
        cpu: cpu
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
    output {
        File paddedCnvBed = "padded_cnv.bed"
    }
}

task SamtoolsDepth {
        input {
            String sampleName
            File alignedBam
            File alignedBai
            File target_bed
            Int minBaseQuality = 20
            Int minMappingQuality = 20
            Int mem_gb = 32
            Int cpu = 4
            Int disk_size_gb = 500
            Boolean use_ssd = true
            String samtools_docker = "euformatics/samtools:1.20"
    }
    command <<<
        # Create directories for input & output
        mkdir input
        mkdir output
        readlink -f ~{alignedBam} > input/bam_path.txt

        # Run samtools depth
        # Counting fragments instead of reads using -s option
        samtools depth \
        -b ~{target_bed} \
        -f input/bam_path.txt \
        --min-BQ ~{minBaseQuality} \
        --min-MQ ~{minMappingQuality} \
        -s \
        -o output/~{sampleName}_samtools.depth

    >>>
    output {
        File depth_profile = "output/~{sampleName}_samtools.depth"
    }
    runtime {
        memory: mem_gb * 1000 + " MB"
        cpu: cpu
        docker: samtools_docker
        disks: "local-disk " + disk_size_gb + if use_ssd then " SSD" else " HDD"
        preemptible: 0
        maxRetries: 3
    }
}

task cnvDepthProfiler{
    input {
            String sampleName
            String cnvProfiler_Docker
            File depthProfile
            File cnvBedFile
            File PaddedcnvBedFile
            Int intervalPadding = 0
            Int mem_gb = 64
            Int cpu = 8
            Int preemptible = 0
            Int disk_size_gb = 500
            Int maxRetries = 1
            Boolean use_ssd = true
        }
        command <<<
            set -e
            mkdir output

            # Run the coverage profile visualization script
            conda run --no-capture-output \
            -n env_viz \
            python3 /BaseImage/CovProfileViz/scripts/CNV_Depth_Profiler.py \
            -c ~{depthProfile} \
            -b ~{cnvBedFile} \
            -n ~{sampleName} \
            -pb ~{PaddedcnvBedFile} \
            -p ~{intervalPadding} \
            -o output

        >>>
        output {
            Array[File] cnv_depth_profile = glob("output/*png")
        }
        runtime {
            memory: mem_gb + " GB"
            cpu: cpu
            docker: cnvProfiler_Docker
            disks: "local-disk " + disk_size_gb + if use_ssd then " SSD" else " HDD"
            preemptible: preemptible
            maxRetries: maxRetries
        }
}

task HeterozygosityCheck{
    input {
        String sampleName
        String cnvProfiler_Docker
        File HG1_vcf_path = "gs://dragenv4_2_validation/dragen_4_2_4/hg19/NA12878/smoke.hard-filtered.vcf.gz"
        File HG2_vcf_path = "gs://dragenv4_2_validation/dragen_4_2_4/hg19/NA24385/NA24385.hard-filtered.vcf.gz"
        File? hardFilteredVcfFile
        File cnvBedFile
        Int mem_gb = 64
        Int cpu = 8
        Int preemptible = 0
        Int disk_size_gb = 500
        Int maxRetries = 1
        Boolean use_ssd = true
    }
    command <<<
        set -e
        mkdir output

        # Run the coverage profile visualization script
        conda run --no-capture-output \
        -n env_viz \
        python3 /BaseImage/CovProfileViz/scripts/CNV_SNP_HET_Profiler.py \
        -v1 ~{hardFilteredVcfFile} \
        -v2 ~{HG1_vcf_path} \
        -v3 ~{HG2_vcf_path} \
        -b ~{cnvBedFile} \
        -n1 ~{sampleName} \
        -n2 HG001 \
        -n3 HG002  \
        -o output

    >>>
    output {
        Array[File] heterozygosity_plot = glob("output/*png")
    }
    runtime {
        memory: mem_gb + " GB"
        cpu: cpu
        docker: cnvProfiler_Docker
        disks: "local-disk " + disk_size_gb + if use_ssd then " SSD" else " HDD"
        preemptible: preemptible
        maxRetries: maxRetries
    }
}







