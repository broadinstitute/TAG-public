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
                referenceFasta = referenceFasta,
                referenceDict = referenceDict,
                referenceFai = referenceFai,
                target_bed = IntervalListToBed.bed_intervals,
                minBaseQuality = MinBaseQuality,
                minMappingQuality = MinMappingQuality
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
            File referenceFasta
            File referenceDict
            File referenceFai
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