version 1.0

workflow CrosscheckFingerprintsWorkflow {
    input {
        File input_bamA
        File input_bamB
        File hapmap
        File reference_fasta
        String output_name
        String picard_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
    }

    call CrosscheckFingerprints {
        input:
            input_bamA = input_bamA,
            input_bamB = input_bamB,
            hapmap = hapmap,
            reference = reference_fasta,
            output_name = output_name,
            picard_docker = picard_docker
    }

    output {
        File crosscheckmetrics = CrosscheckFingerprints.crosscheckmetrics
    }

    meta {
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        description: "CrossCheck Fingerprints WDL"
    }
}

task CrosscheckFingerprints {
    input {
        File input_bamA
        File input_bamB
        File hapmap
        File reference
        String output_name
        String picard_docker
        Int disk_gb = 200
        Int memory_gb = 16
    }

    command <<<
        set -e
        mkdir temp_dir
        java -jar /usr/gitc/picard.jar CrosscheckFingerprints \
            -I ~{input_bamA} \
            -I ~{input_bamB} \
            -HAPLOTYPE_MAP ~{hapmap} \
            -CROSSCHECK_BY SAMPLE \
            -EXPECT_ALL_GROUPS_TO_MATCH true \
            -O ~{output_name}.crosscheck_metrics \
            -REFERENCE_SEQUENCE ~{reference} \
            -TMP_DIR ./temp_dir
    >>>

    runtime {
        docker: picard_docker
        disks: "local-disk " + disk_gb + " HDD"
        memory: memory_gb + "GB"
        cpu: "1"
    }

    output {
        File crosscheckmetrics = "~{output_name}.crosscheck_metrics"
    }
}

