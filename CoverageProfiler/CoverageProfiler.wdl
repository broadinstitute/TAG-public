version 1.0

workflow coverageProfile {
    input {
        String sampleName
        String coverageTool = "Samtools"
        File alignedBam
        File alignedBamIndex
        File referenceFasta
        File referenceDict
        File referenceFai
        File intervals
        File? interval_GCcontent_track
        Int MinBaseQuality = 20
        Int MinMappingQuality = 20
        Boolean visualise_coverage = false
    }
    if (coverageTool =="Samtools") {
        call IntervalListToBed {
            input:
                intervals = intervals
        }
        call SamtoolsDepth {
            input:
                sampleName = sampleName,
                alignedBam = alignedBam,
                alignedBamIndex = alignedBamIndex,
                target_bed = IntervalListToBed.bed_intervals,
                minBaseQuality = MinBaseQuality,
                minMappingQuality = MinMappingQuality
        }
        if (visualise_coverage) {
            call CovProfileViz {
                input:
                    sampleName = sampleName,
                    SamtoolsDepthProfile = SamtoolsDepth.depth_profile,
                    GCcontentTrack = interval_GCcontent_track
            }
            }

        }
    if (coverageTool == "DepthOfCoverage") {
        call DepthOfCoverage {
            input:
                sampleName = sampleName,
                alignedBam = alignedBam,
                alignedBamIndex = alignedBamIndex,
                referenceFasta = referenceFasta,
                referenceDict = referenceDict,
                referenceFai = referenceFai,
                intervals = intervals,
                minBaseQuality = MinBaseQuality,
                minMappingQuality = MinMappingQuality
        }
    }
    output {
        File? DepthOfCoverageIntervalCov = DepthOfCoverage.sample_interval_summary
        Float? DepthOfCoverageMeanCoverage = DepthOfCoverage.mean_coverage
        File? SamtoolsDepthProfile = SamtoolsDepth.depth_profile
        File? SamtoolsCovProfilePlot = CovProfileViz.cov_profile_plot
        Float? SamtoolsAvgCovMean = CovProfileViz.avg_cov_mean
    }
    meta {
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        description:"Calculates the depth of coverage of an input sample and visualize it. Currently supports two tools: Samtools and DepthOfCoverage. The visualisation is only available for Samtools and Exome."

    }
}
    task DepthOfCoverage {
        input {
            String sampleName
            File alignedBam
            File alignedBamIndex
            File referenceFasta
            File referenceDict
            File referenceFai
            File intervals
            Int minBaseQuality
            Int minMappingQuality
            Int mem_gb = 32
            Int cpu = 4
            String gatk_docker = "broadinstitute/gatk:4.5.0.0"
        }
            Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
            Int command_mem_mb = machine_mem_mb - 1000
        command <<<
            # Create directories for output
            mkdir output

            # Run DepthOfCoverage
            # Count all reads independently (even if from the same fragment). The only option supported by GATK 4.5.0.0.
            gatk --java-options "-Xmx~{command_mem_mb}m" DepthOfCoverage \
                -L ~{intervals} \
                --input ~{alignedBam} \
                --read-index ~{alignedBamIndex} \
                --reference ~{referenceFasta} \
                --minimum-mapping-quality ~{minMappingQuality} \
                --min-base-quality ~{minBaseQuality} \
                --count-type COUNT_READS \
                --output output/~{sampleName}

            cat output/~{sampleName}.sample_interval_summary | awk 'BEGIN {FS = ","}{print $3}' | tail -n 1 > output/mean_coverage.txt
        >>>
        output {
            File sample_interval_summary = "output/~{sampleName}.sample_interval_summary"
            Float mean_coverage = read_float("output/mean_coverage.txt")
        }
        runtime {
            memory: machine_mem_mb + " MB"
            cpu: select_first([cpu, 1])
            docker: gatk_docker
            disks: "local-disk 500 SSD"
        }
}
    task IntervalListToBed {
        input {
            File intervals
        }
        command <<<
            # Create directories for output
            mkdir output

            # Convert interval list to BED file
            # This is necessary because samtools depth requires a BED file as input
            gatk IntervalListToBed \
                --INPUT ~{intervals} \
                --OUTPUT output/intervals.bed
        >>>
        output {
            File bed_intervals = "output/intervals.bed"
        }
        runtime {
            docker: "broadinstitute/gatk:4.5.0.0"
            memory: "1 GB"
            cpu: 1
        }
    }

    task SamtoolsDepth {
        input {
            String sampleName
            File alignedBam
            File alignedBamIndex
            File target_bed
            Int minBaseQuality
            Int minMappingQuality
            Int mem_gb = 32
            Int cpu = 4
            Int disk_size_gb = 500
            Int preemptible = 1
            Int maxRetries = 3
            String samtools_docker = "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
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
        memory: select_first([mem_gb, 7]) * 1000 + " MB"
        cpu: select_first([cpu, 1])
        docker: samtools_docker
        disks: "local-disk ~{disk_size_gb} SSD"
        preemptible: preemptible
        maxRetries: maxRetries
    }
}
    task CovProfileViz {
        input {
            File SamtoolsDepthProfile
            File? GCcontentTrack
            String sampleName
            String CovProfileViz_docker = "us-central1-docker.pkg.dev/tag-team-160914/gptag-dockers/covprofileviz:0.0.0"
            Int mem_gb = 32
            Int cpu = 4
            Int preemptible = 1
            Int MaxRetries = 3
            Int disk_size_gb = 500
        }
        command <<<
            set -e
            mkdir output
            # Run the coverage profile visualization script
            conda run --no-capture-output \
            -n env_viz \
            python3 /BaseImage/CovProfileViz/scripts/plot_samtoolsDepths_by_chr.py \
            -s ~{sampleName} \
            -d ~{SamtoolsDepthProfile} \
            -g ~{GCcontentTrack} \
            -o output

            mv output/*_samtools_cov_with_gc.png output/~{sampleName}Sample_Cov_profile.png
            mv output/*_avg_cov_mean.txt output/~{sampleName}_Avg_Cov_mean.txt
        >>>
        output {
            File cov_profile_plot = "output/~{sampleName}Sample_Cov_profile.png"
            Float avg_cov_mean = read_float("output/~{sampleName}_Avg_Cov_mean.txt")
        }
        runtime {
            memory: mem_gb + " GB"
            cpu: select_first([cpu, 1])
            docker: CovProfileViz_docker
            disks: "local-disk ~{disk_size_gb} SSD"
            preemptible: preemptible
            maxRetries: MaxRetries
        }
    }