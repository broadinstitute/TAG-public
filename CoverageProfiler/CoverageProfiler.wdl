version 1.0

workflow coverageProfile {
    input {
        String sampleName
        File alignedBam
        File referenceFasta
        File referenceDict
        File referenceFai
        File intervals
    }
    call DepthOfCoverage {
        input:
            sampleName = sampleName,
            alignedBam = alignedBam,
            referenceFasta = referenceFasta,
            referenceDict = referenceDict,
            referenceFai = referenceFai,
            intervals = intervals
    }
    output {
        File coveragebyInterval = DepthOfCoverage.sample_interval_summary
        Float meanCoverage = DepthOfCoverage.mean_coverage
    }
    meta {
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        description: "Calculates the depth of coverage of an input sample using GATK's DepthOfCoverage tool."

    }
}
    task DepthOfCoverage {
    input {
        String sampleName
        File alignedBam
        File referenceFasta
        File referenceDict
        File referenceFai
        File intervals
        Int? mem_gb
        Int? cpu
        String gatk_docker = "broadinstitute/gatk:4.5.0.0"
    }
        Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
        Int command_mem_mb = machine_mem_mb - 1000
    command <<<
        # Create directories for input, output, reference localization
        mkdir input
        mkdir output
        mkdir reference

        # Localize Reference Files and BAM
        mv ~{referenceFasta} reference/reference.fasta
        mv ~{referenceDict} reference/reference.dict
        mv ~{referenceFai} reference/reference.fai
        mv ~{alignedBam} input/~{sampleName}.bam

        # Index BAM file
        gatk BuildBamIndex \
            --INPUT input/~{sampleName}.bam

        # Run DepthOfCoverage
        gatk --java-options "-Xmx~{command_mem_mb}m" DepthOfCoverage \
            -L ~{intervals} \
            --input input/~{sampleName}.bam \
            --reference reference/reference.fasta \
            --output output/~{sampleName}

        cat output/~{sampleName}.sample_interval_summary | awk 'BEGIN {FS = ","}{print $3}' | tail -n 1 > output/mean_coverage.txt
    >>>
    output {
        File sample_interval_summary = "output/${sampleName}.sample_interval_summary"
        Float mean_coverage = read_float("output/mean_coverage.txt")
    }
    runtime {
        memory: machine_mem_mb + " MB"
        cpu: select_first([cpu, 1])
        docker: gatk_docker
        disks: "local-disk 500 HDD"
    }
}