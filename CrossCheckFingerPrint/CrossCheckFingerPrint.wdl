version 1.0

workflow CrosscheckFingerprintsWorkflow {
    input {
        Array[File] input_bams
        Boolean fail_on_mismatch = false
        Int mismatch_rc = if fail_on_mismatch then 1 else 2
        File hapmap
        File reference_fasta
        String output_name
        String picard_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
    }

    call CrosscheckFingerprints {
        input:
            input_bams = input_bams,
            hapmap = hapmap,
            reference = reference_fasta,
            output_name = output_name,
            picard_docker = picard_docker,
            mismatch_rc = mismatch_rc
    }

    output {
        File crosscheckmetrics = CrosscheckFingerprints.crosscheckmetrics
        String crosscheck_results = CrosscheckFingerprints.crosscheck_results
    }

    meta {
        author: "Yueyao Gao, Micah Rickles-Young"
        email: "tag@broadinstitute.org"
        description: "CrossCheck Fingerprints WDL"
    }
}

task CrosscheckFingerprints {
    input {
        Array[File] input_bams
        File hapmap
        File reference
        String output_name
        String picard_docker
        Int mismatch_rc
        String extra_args = ""
        Int disk_gb = 200
        Int memory_gb = 16
    }

    command <<<
        set -e
        mkdir temp_dir
        java -jar /usr/gitc/picard.jar CrosscheckFingerprints \
        I=~{sep=' I=' input_bams} \
        HAPLOTYPE_MAP=~{hapmap} \
        CROSSCHECK_BY=SAMPLE \
        O=~{output_name}.crosscheck_metrics.txt \
        REFERENCE_SEQUENCE=~{reference} \
        EXIT_CODE_WHEN_MISMATCH=~{mismatch_rc} \
        ~{extra_args} \
        TMP_DIR=./temp_dir
    >>>

    runtime {
        docker: picard_docker
        disks: "local-disk " + disk_gb + " HDD"
        continueOnReturnCode: [0, 2]
        memory: memory_gb + "GB"
        cpu: "1"
    }

    output {
        File crosscheckmetrics = "~{output_name}.crosscheck_metrics.txt"
        String crosscheck_results = if read_int("rc") == 0 then "PASS" else "FAIL"
    }
}
