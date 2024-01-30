version 1.0


workflow demux_CODEC {
    input {
        Array[File] fastq1
        Array[File] fastq2
        Array[String] batch_ids
        Int num_parallel
        Array[File] sample_sheet
        String workspace
    }
    scatter (batch_index in range(length(batch_ids))) {
        call SplitFastq1 {
            input:
                batch_id =batch_ids[batch_index],
                fastq_read1 = fastq1[batch_index],
                nsplit = num_parallel
        }
        call SplitFastq2 {
            input:
                batch_id =batch_ids[batch_index],
                fastq_read2 = fastq2[batch_index],
                nsplit = num_parallel
        }
        scatter (split_index in range(num_parallel)) {
            call Demux {
                input:
                    paired_fastqs = (SplitFastq1.split_read1[split_index], SplitFastq2.split_read2[split_index]),
                    sample_sheet = sample_sheet[batch_index],
                    batch_id = batch_ids[batch_index],
                    split = split_index
        }
    }
    }
    Array[File] demultiplexed_fastq1_all_lanes = flatten(flatten(Demux.demux_read1))
    Array[File] demultiplexed_fastq2_all_lanes = flatten(flatten(Demux.demux_read2))
    
    call GroupFastqFiles as GroupFastq1 {
      input:
      fastq_files = demultiplexed_fastq1_all_lanes
    }
    call GroupFastqFiles as GroupFastq2 {
      input:
      fastq_files = demultiplexed_fastq2_all_lanes
    }
    scatter (txt_file in GroupFastq1.grouped_fastq_lists) {
        call concatenate_and_compress_fastq as concatenate_and_compress_fastq1 {
            input:
                grouped_fastq_list = txt_file,
                output_file = basename(txt_file, ".txt") + ".1.gz"
        }
    }

    scatter (txt_file in GroupFastq2.grouped_fastq_lists) {
        call concatenate_and_compress_fastq as concatenate_and_compress_fastq2 {
            input:
                grouped_fastq_list = txt_file,
                output_file = basename(txt_file, ".txt") + ".2.gz"
        }
    }
    
    call UploadDataTable {
        input:
           fastq1_files = concatenate_and_compress_fastq1.output_fastq_gz,
           fastq2_files = concatenate_and_compress_fastq2.output_fastq_gz,
           workspace = workspace

    } 

    output {
        Array[File] merged_fastq1_per_sample = concatenate_and_compress_fastq1.output_fastq_gz
        Array[File] merged_fastq2_per_sample = concatenate_and_compress_fastq2.output_fastq_gz
        File data_table = UploadDataTable.data_table
    }
}


task SplitFastq1 {
    input {
        File fastq_read1
        Int nsplit
        String batch_id
        Int memory = 128
        Int disk_size = 1024
 
    }

    command <<<
        set -e
        
        zcat ~{fastq_read1} | /CODECsuite/snakemake/script/fastqsplit.pl ~{batch_id}_split_r1 ~{nsplit}

    >>>

    output {
        Array[File] split_read1 = glob("~{batch_id}_split_r1.*.fastq")
    }

    runtime {
        docker: "us.gcr.io/tag-team-160914/codec:v1"
        memory: memory + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task SplitFastq2 {
    input {
        File fastq_read2
        Int nsplit
        String batch_id
        Int memory = 128
        Int disk_size = 1024
    }

    command <<<
        set -e
        zcat ~{fastq_read2} | /CODECsuite/snakemake/script/fastqsplit.pl ~{batch_id}_split_r2 ~{nsplit} 
        
    >>>

    output {
        Array[File] split_read2 = glob("~{batch_id}_split_r2.*.fastq")
    }

    runtime {
        docker: "us.gcr.io/tag-team-160914/codec:v1"
        memory: memory + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task Demux {
    input {
        Pair[File, File] paired_fastqs
        File sample_sheet
        Int memory = 64
        Int disk_size = 100
        Int split
        String batch_id
    }
    
    command {
        set -e 
        
        /CODECsuite/build/codec demux -1 ${paired_fastqs.left} -2 ${paired_fastqs.right} -p ${sample_sheet} -o ${batch_id}_split.${split} > ${batch_id}_split.${split}.log
    }
    runtime {
        docker: "us.gcr.io/tag-team-160914/codec:v1"
        memory: memory + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
    output {
        Array[File] demux_read1 = glob("~{batch_id}_split.~{split}.*.1.fastq.gz")
        Array[File] demux_read2 = glob("~{batch_id}_split.~{split}.*.2.fastq.gz")
        File log = "~{batch_id}_split.~{split}.log"
    }
}
task GroupFastqFiles {
  input {
    Array[File] fastq_files
    Int memory = 64
    Int disk_size = 100
  }

  command <<<
    set -e

    for fastq_file in ~{sep=' ' fastq_files}; do
        sample_name=$(basename "$fastq_file" | cut -d '.' -f 3)
        echo "$fastq_file" >> "${sample_name}.txt"
    done
    
    for txt_file in *.txt; do
    sed -i 's|/cromwell_root/|gs://|g' "$txt_file"
    done
  >>>

  output {
    Array[File] grouped_fastq_lists = glob("*.txt")
  }

  runtime {
    docker: "us.gcr.io/tag-team-160914/codec:v1"
    disks: "local-disk " + disk_size + " HDD"
    memory: memory + " GB"
  }
}

task concatenate_and_compress_fastq {
    input {
        File grouped_fastq_list
        String output_file
        Int memory = 64
        Int disk_size = 64
    }

    command {
        set -e
        while IFS= read -r line; do
        gsutil cp "$line" temp.fastq.gz
        
        # Check if the file was successfully copied
        if [ -f temp.fastq.gz ]; then
            gunzip -c temp.fastq.gz
        else
            echo "Failed to copy file from $line"
            exit 1
        fi
    done < ~{grouped_fastq_list} | gzip > ~{output_file}
    }
    output {
        File output_fastq_gz = output_file
    }

    runtime {
        docker: "us.gcr.io/tag-team-160914/gcloudsdk"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory + " GB"
    }
}


task UploadDataTable {
    input {
        Array[String] fastq1_files
        Array[String] fastq2_files
        String output_data_table = "sample.tsv"
        String workspace
        String namespace = "broadtagteam"
        Int memory = 32
        Int disk_size = 32
    }

    command <<<
        python3 <<CODE
        import os
        import firecloud.api as fapi
        import argparse

        fastq1_files = "~{sep=' ' fastq1_files}".split()
        fastq2_files = "~{sep=' ' fastq2_files}".split()
        output_data_table = "~{output_data_table}"

        if len(fastq1_files) != len(fastq2_files):
            raise ValueError("The length of fastq1_files and fastq2_files must be the same.")

        with open(output_data_table, "w") as out_file:
            # Write the header
            out_file.write("entity:sample_id\tfastq1\tfastq2\n")

            for fq1, fq2 in zip(fastq1_files, fastq2_files):
                sample_id = os.path.basename(fq1).split('.')[0]
                out_file.write(f"{sample_id}\t{fq1}\t{fq2}\n")

        upload_sample_response = fapi.upload_entities_tsv( "~{namespace}", "~{workspace}",
                                                        output_data_table, model='flexible')
        if "message" not in upload_sample_response:
            print("Sample table was uploaded successfully!")
        else:
            print("Error uploading sample table:",
                upload_sample_response["message"])

        CODE
    >>>

    output {
        File data_table = output_data_table
        String upload_stats = read_string(stdout())
    }

    runtime {
        docker: "us.gcr.io/tag-team-160914/metadata_upload"
        disks: "local-disk " + disk_size + " HDD"
        memory: memory + " GB"
    }
}