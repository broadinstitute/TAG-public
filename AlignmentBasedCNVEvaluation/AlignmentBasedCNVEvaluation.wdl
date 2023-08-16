version 1.0

    workflow AlignmentBasedCNVEval {
        input {
            File cnv_vcf
            File reference_fasta
            File reference_last_database
            File T2T_last_database
            String last_docker
            String analysis_docker
        }
        call extract_cnv {
            input:
                cnv_vcf = cnv_vcf
        }
        scatter (cnv_event_chunk in extract_cnv.cnv_intervals_chunks) {
            call get_sequence_from_interval{
                input:
                    cnv_event_chunk = cnv_event_chunk,
                    reference_fasta = reference_fasta,
                    analysis_docker = analysis_docker
            }
            Array[File] interval_sequence_fasta = get_sequence_from_interval.interval_sequence_fasta

            call last_align {
                input:
                    interval_sequence_fasta = interval_sequence_fasta,
                    reference_last_database = reference_last_database,
                    T2T_last_database = T2T_last_database,
                    last_docker = last_docker
            }
            Array[File] ref_lastal_alignment = last_align.ref_lastal_alignment
            Array[File] t2t_lastal_alignment = last_align.t2t_lastal_alignment

            call evaluate_cnv {
                input:
                    ref_lastal_alignment = ref_lastal_alignment,
                    t2t_lastal_alignment = t2t_lastal_alignment,
                    analysis_docker = analysis_docker,
                    uncompressed_cnv_vcf = extract_cnv.uncompressed_cnv_vcf,

            }
        }
        # Flatten the array of arrays
        Array[File] cnv_copy_number = flatten(evaluate_cnv.cnv_copy_number)
        call gather_results {
            input:
                analysis_docker = analysis_docker,
                cnv_copy_number = cnv_copy_number
        }
        output {
            File cnv_eval_alignment_output = gather_results.gathered_alignment_output
        }
    }

    task extract_cnv {
            input {
                File cnv_vcf
                Int cpu_num = 1
                Int mem_size = 4
                Int disk_size = 10
            }
            command <<<
                set -e

                filename=$(basename ~{cnv_vcf})
                # Unzip the input file if it is gzipped
                if test "${filename##*.}" = "gz"; then
                  gunzip -c ~{cnv_vcf} > temp.vcf
                else
                  cp ~{cnv_vcf} temp.vcf
                fi

                # Extract the PASS events from the DRAGEN CNV VCF and save them to a file
                cat temp.vcf | grep -v '#' | awk '$7 == "PASS" {print}' | awk '{print $3}' | awk 'BEGIN {FS=":"} {print $1":"$2"\t"$3":"$4}' > pass_cnv_event.txt

                # Split CNV intervals to individual files
                split -l 100 -d --additional-suffix=_interval.txt pass_cnv_event.txt pass_cnv_chunk_

                # Print the number of GAINS and LOSSES that DRAGEN found
                cat temp.vcf | grep -v '#' | awk '$7 == "PASS" {print}' | awk '{print $3}' | awk 'BEGIN {FS=":"} {print $2}' | sort | uniq -c | sort
            >>>
            output {
                File uncompressed_cnv_vcf = "temp.vcf"
                File cnv_events_file = "pass_cnv_event.txt"
                Array[File] cnv_intervals_chunks = glob("pass_cnv_chunk_*")

            }
            runtime {
                docker: "us.gcr.io/broad-dsde-methods/liquidbiopsy:0.0.3.7"
                bootDiskSizeGb: 12
                cpu: cpu_num
                memory: mem_size + " GB"
                disks: "local-disk " + disk_size + " HDD"
                preemptible: 2
                #maxRetries: 3
            }
        }

    task get_sequence_from_interval {
        input {
            File cnv_event_chunk
            File reference_fasta
            String analysis_docker
        }
        command <<<
            set -e

            # Obtain sequence fasta file for each CNV interval
            for i in `cat ~{cnv_event_chunk} | awk '{print $2}'`; do python3 /blastn/getSeq.py -i ${i} -r ~{reference_fasta}; done
        >>>
        output {
            Array[File] interval_sequence_fasta = glob("*_seq.fasta")
        }

        runtime {
                docker: analysis_docker
                bootDiskSizeGb: 12
                cpu: 1
                memory: "4 GB"
                disks: "local-disk 100 HDD"
                preemptible: 2
                maxRetries: 3
        }
    }

    task last_align {
        input {
            Array [File] interval_sequence_fasta
            File reference_last_database
            File T2T_last_database
            String last_docker
            Int diskSizeGB = 1000
        }
        command <<<
            set -e

            # Import the resource monitor script
            export TMPDIR=/tmp
            wget -O monitor_script.sh https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh
            chmod a+rx monitor_script.sh
            ./monitor_script.sh > monitoring.log &

            # Extract last database from tar files
            mkdir -p /cromwell_root/reference_database
            mkdir -p /cromwell_root/t2t_database

            tar -xvf ~{reference_last_database} -C /cromwell_root/reference_database/
            tar -xvf ~{T2T_last_database} -C /cromwell_root/t2t_database/

            # Basename for the blast database
            # If there are dot in the name, it will cause error in extracting the basename
            reference_db_path=$(echo "`readlink -f /cromwell_root/reference_database/*`/`basename /cromwell_root/reference_database/*/*.prj | cut -d '.' -f 1`")
            t2t_db_path=$(echo "`readlink -f /cromwell_root/t2t_database/*`/`basename /cromwell_root/t2t_database/*/*.prj | cut -d '.' -f 1`")

            # Run lastal
            cd /cromwell_root
            for i in `ls ~{sep=" " interval_sequence_fasta}`; do
                interval_name=$(basename ${i} "_seq.fasta")
                /last-1460/bin/lastal5 $reference_db_path ${i} -v -P 0 -l 30 -f BlastTab > /cromwell_root/${interval_name}_ref_lastal_alignment.txt
                /last-1460/bin/lastal5 $t2t_db_path ${i} -v -P 0 -l 30 -f BlastTab > /cromwell_root/${interval_name}_t2t_lastal_alignment.txt
            done

        >>>
        output {
            Array[File] ref_lastal_alignment = glob("*_ref_lastal_alignment.txt")
            Array[File] t2t_lastal_alignment = glob("*_t2t_lastal_alignment.txt")
        }

        runtime {
                docker: last_docker
                bootDiskSizeGb: 12
                cpu: 16
                memory: "128 GB"
                disks: "local-disk ~{diskSizeGB} HDD"
                preemptible: 0
                maxRetries: 0
        }
    }

    task evaluate_cnv {
        input {
            Array[File] ref_lastal_alignment
            Array[File] t2t_lastal_alignment
            File uncompressed_cnv_vcf
            String analysis_docker
        }
        command <<<
            set -e
            # Using bash to process the file arrays and find matching filenames
            for ref_lastal_output in ~{sep=' ' ref_lastal_alignment}; do
                for t2t_lastal_output in ~{sep=' ' t2t_lastal_alignment}; do
                    # Extract the filename without path
                    filename1=$(basename $ref_lastal_output "_ref_lastal_alignment.txt")
                    filename2=$(basename $t2t_lastal_output "_t2t_lastal_alignment.txt")

                    # Check if the filenames match
                    if [ "$filename1" == "$filename2" ]; then
                        python3 /blastn/lastal_parser.py -r $ref_lastal_output -t $t2t_lastal_output
                        cat ~{uncompressed_cnv_vcf} | grep "$filename1" | awk '{print $5}' > cnv_event_type.txt
                        paste cnv_event_type.txt ${filename1}_copy_number.txt > ${filename1}_annoed_copy_number.txt
                    fi
                done
            done
    >>>
        output {
            Array[File] cnv_copy_number = glob("*_annoed_copy_number.txt")
        }
        runtime {
                docker: analysis_docker
                bootDiskSizeGb: 12
                cpu: 1
                memory: "32 GB"
                disks: "local-disk 300 HDD"
                preemptible: 0
                maxRetries: 0
        }
    }

    task gather_results {
        input{
            Array[File] cnv_copy_number
            String analysis_docker
         }
        command <<<
        set -e

        cat ~{sep=" " cnv_copy_number} >> aggregated_cnv_alignment_output.txt
            >>>
        output {
            File gathered_alignment_output = "aggregated_cnv_alignment_output.txt"
        }
        runtime {
                docker: analysis_docker
                bootDiskSizeGb: 12
                cpu: 1
                memory: "32 GB"
                disks: "local-disk 100 SSD"
                preemptible: 2
                maxRetries: 3
        }
    }