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
        File interval_GCcontent_track
        Int MinBaseQuality = 20
        Int MinMappingQuality = 20
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
        call CovProfileViz {
            input:
                sampleName = sampleName,
                SamtoolsDepthProfile = SamtoolsDepth.depth_profile,
                GCcontentTrack = interval_GCcontent_track
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
        Float? SamtoolsAvgChrCovStd = CovProfileViz.avg_chr_cov_std
        File? SamtoolsAvgChrCovPerChr = CovProfileViz.avg_chr_cov_per_chr
        Float? SamtoolsAvgCovMean = CovProfileViz.avg_cov_mean
        File? SamtoolsAvgChrCovPerChrPlot = CovProfileViz.avg_chr_cov_per_chr_plot
    }
    meta {
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        description: "Calculates the depth of coverage of an input sample"

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
            Int? mem_gb
            Int? cpu
            String gatk_docker = "broadinstitute/gatk:4.5.0.0"
        }
            Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
            Int command_mem_mb = machine_mem_mb - 1000
        command <<<
            # Create directories for output
            mkdir output

            # Run DepthOfCoverage
            gatk --java-options "-Xmx~{command_mem_mb}m" DepthOfCoverage \
                -L ~{intervals} \
                --input ~{alignedBam} \
                --read-index ~{alignedBamIndex} \
                --reference ~{referenceFasta} \
                --minimum-mapping-quality ~{minMappingQuality} \
                --min-base-quality ~{minBaseQuality} \
                --count-type COUNT_READS \ # Count all reads independently (even if from the same fragment). The only option supported by GATK 4.5.0.0.
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
            Int? mem_gb
            Int? cpu
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
        memory: select_first([mem_gb, 7]) * 1000 + " MB"
        cpu: select_first([cpu, 1])
        docker: samtools_docker
        disks: "local-disk 500 SSD"
    }
}
    task CovProfileViz {
        input {
            File SamtoolsDepthProfile
            File GCcontentTrack
            String sampleName
            String CovProfileViz_docker = "us-central1-docker.pkg.dev/tag-team-160914/gptag-dockers/covprofileviz:0.0.0"
            Int mem_gb = 32
            Int? cpu
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

            mv output/*_samtools_cov_with_gc.png output/sample_coverage_profile.png
            mv output/*_avg_cov_per_chr.png output/avg_cov_per_chr.png
            mv output/*_avg_cov_std.txt output/per_chr_cov_std.txt
            mv output/*_avg_cov_per_chr.csv output/per_chr_avg_cov.csv
            mv output/*_avg_cov_mean.txt output/avg_cov_mean.txt
        >>>
        output {
            File cov_profile_plot = "output/sample_coverage_profile.png"
            File avg_chr_cov_per_chr_plot = "output/avg_cov_per_chr.png"
            Float avg_chr_cov_std = read_float("output/per_chr_cov_std.txt")
            File avg_chr_cov_per_chr = "output/per_chr_avg_cov.csv"
            Float avg_cov_mean = read_float("output/avg_cov_mean.txt")
        }
        runtime {
            memory: select_first([mem_gb, 7]) * 1000 + " MB"
            cpu: select_first([cpu, 1])
            docker: CovProfileViz_docker
        }
    }