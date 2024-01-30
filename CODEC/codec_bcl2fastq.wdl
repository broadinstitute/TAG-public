version 1.0

workflow codec_bcl2fastq {
    input {
        String input_bcl_directory
        String output_directory
        Int read_threads 
        Int write_threads
        String memory = "128G"
        Int disk_space = 2000
        Int preemptible = 3
        Boolean delete_input_bcl_directory
    }

    call run_bcl2fastq {
        input:
            input_bcl_directory = sub(input_bcl_directory, "/+$", ""),
            output_directory = sub(output_directory, "/+$", ""),
            read_threads = read_threads,
            write_threads = write_threads,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible,  
            delete_input_bcl_directory = delete_input_bcl_directory
    }

    output {
        String job_stat = run_bcl2fastq.job_stat
    }
}

task run_bcl2fastq {
    input {
        String input_bcl_directory
        String output_directory
        String memory
        Int disk_space
        Int preemptible
        Int read_threads = 4
        Int write_threads = 4
        String run_id = basename(input_bcl_directory)
        Boolean delete_input_bcl_directory
        Int num_cpu = 2
    }

    command {
        set -e
        export TMPDIR=/tmp

        gsutil -q -m cp -r ${input_bcl_directory} .

        cd ${run_id}

        bcl2fastq \
        --output-dir out \
        -w ~{write_threads} \
        -r ~{read_threads}

        cd out

        gsutil -q -m cp -r . ${output_directory}/${run_id}_fastqs/

        if [ $? -eq 0 ]; then
            echo "Fastq files generated and copied successfully to ${output_directory}/${run_id}_fastqs/."
        else
            echo "Error in copying files."
        fi

        if [ "${delete_input_bcl_directory}" = "true" ]; then
            gsutil -q -m rm -r "${input_bcl_directory}"
            if [ $? -eq 0 ]; then
                echo "Input bcl folder has been successfully deleted!"
            else
                echo "Unable to delete BCL directory"
            fi
        fi
    }

    output {
        String job_stat = read_string(stdout())
    }

 	runtime {
		docker: "us.gcr.io/tag-team-160914/bcl2fastq_codec:v1"
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
}
}