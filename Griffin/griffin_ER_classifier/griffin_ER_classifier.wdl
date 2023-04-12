version 1.0

workflow ER_status_classification{
    input {
        String sample_name
        File samples_yaml
        File GC_corrected_cov
        File ER_status_scaler
        File ER_status_model
        String griffin_docker

    }
    call ER_classifier {
        input:
            sample_name = sample_name,
            samples_yaml = samples_yaml,
            GC_corrected_cov = GC_corrected_cov,
            ER_status_scaler = ER_status_scaler,
            ER_status_model = ER_status_model,
            griffin_docker = griffin_docker
    }
    output {
        Int ER_status = ER_classifier.ER_status
        Float ER_prob = ER_classifier.ER_prob
    }
    meta {
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        description: "Griffin: ER subtyping for MBC samples"
    }
}

task ER_classifier {
    input {
        String sample_name
        File samples_yaml
        File GC_corrected_cov
        File ER_status_scaler
        File ER_status_model
        String griffin_docker
        Int? cpu
        Int? mem
        Int? disk_space
        # If cpu, mem, and disk size were not specified, use 8 cores, 10GB, and 100 GB as default
        Int cpu_num = select_first([cpu, 1])
        Int mem_size = select_first([mem, 4])
        Int disk_size = select_first([disk_space,100])
    }
    command <<<
        set -e
        # Run griffin_classifer_script to calculate ER status
        conda run --no-capture-output \
        -n griffin_env \
        python3 /BaseImage/Griffin/scripts/griffin_ER_classifier.py \
        --samples_name ~{sample_name} \
        --samples_yaml ~{samples_yaml} \
        --gc_corrected_coverage ~{GC_corrected_cov} \
        --scaler ~{ER_status_scaler} \
        --model ~{ER_status_model} \

        awk -F $'\t' '{print $2}' ~{sample_name}_ER_status_prediction.tsv | sed -n 2p > ER_status.txt
        awk -F $'\t' '{print $3}' ~{sample_name}_ER_status_prediction.tsv | sed -n 2p > ER_prob.txt
    >>>
    runtime {
        docker: griffin_docker
        bootDiskSizeGb: 12
        cpu: cpu_num
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        maxRetries: 3
        }
    output {
        Int ER_status = read_int("ER_status.txt")
        Float ER_prob = read_float("ER_prob.txt")
    }
}
