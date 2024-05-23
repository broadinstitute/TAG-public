version 1.0

workflow codec_bcl2fastq {
    input {
        String input_bcl_directory
        String output_directory
        Int read_threads 
        Int write_threads
        String memory = "128G"
        Int preemptible = 3
        Boolean delete_input_bcl_directory
    }
    call GetBclBucketSize {
        input:
            input_bcl_directory = input_bcl_directory,
            memory = 16
    }

    call run_bcl2fastq {
        input:
            input_bcl_directory = sub(input_bcl_directory, "/+$", ""),
            output_directory = sub(output_directory, "/+$", ""),
            read_threads = read_threads,
            write_threads = write_threads,
            memory = memory,
            input_bcl_size = GetBclBucketSize.input_bcl_size,
            preemptible = preemptible,  
            delete_input_bcl_directory = delete_input_bcl_directory
    }

    output {
        String job_stat = run_bcl2fastq.job_stat
    }
}

task GetBclBucketSize {
  input {
    String input_bcl_directory
    Int memory = 16
  }

  command <<< 
    gsutil du -s ${input_bcl_directory} | awk '{print $1/1024/1024/1024}' > input_bcl_size.txt
  >>>

  output {
    Float input_bcl_size = read_float("input_bcl_size.txt")
  }

  runtime {
    docker: "us.gcr.io/tag-public/bcl2fastq_codec:v1"
    memory: memory + " GB"
	disks: "local-disk 16 HDD"
  }
}

task run_bcl2fastq {
    input {
        String input_bcl_directory
        String output_directory
        String memory
        Float input_bcl_size 
        Int disk_space = ceil(input_bcl_size * 1.5) + select_first([extra_disk, 0])
        Int? extra_disk
        Int read_threads = 4
        Int write_threads = 4
        String run_id = basename(input_bcl_directory)
        Boolean delete_input_bcl_directory
        Int num_cpu = 2
        Int preemptible
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
		docker: "us.gcr.io/tag-public/bcl2fastq_codec:v1"
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
}
}