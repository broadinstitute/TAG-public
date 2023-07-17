version 1.0

    workflow cnv_blast{
        input {
            File cnv_vcf
            File reference_fasta
            File reference_blast_database
            File T2T_blast_database
            String blast_docker = "us.gcr.io/tag-team-160914/cnv_blastn:1.0.1"
        }
        call extract_cnv {
            input:
                cnv_vcf = cnv_vcf
        }

        scatter (cnv_event_file in extract_cnv.cnv_intervals_files) {
            call length_check {
                input:
                    cnv_interval = cnv_event_file
            }

            call blastn {
                input:
                    docker = blast_docker,
                    cnv_length = length_check.cnv_event_length,
                    cnv_interval = cnv_event_file,
                    cnv_vcf = cnv_vcf,
                    reference_fasta = reference_fasta,
                    reference_blast_database = reference_blast_database,
                    T2T_blast_database = T2T_blast_database
            }
        }
        Array[File] blasted_interval_results = blastn.blasted_cnv_interval
        call gather_results {
            input:
                cnv_blasted_results = blasted_interval_results
        }

        output {
            File cnv_events_file = extract_cnv.cnv_events_file
            Array[File] cnv_events = extract_cnv.cnv_intervals_files
            File cnv_intervals_file = gather_results.gathered_intervals_txt
        }
        meta {
            author: "Yueyao Gao"
            email: "tag@broadinstitute.org"
            description: "This WDL script validates DRAGEN CNV calls by cross-referencing them against a ground truth, and utilizes BLAST to search for the CNV sequence in the HG002 T2T and hg38 references, outputting the number of significant matches in each."
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

                # Save the CNV intervals to individual files
                for i in `cat pass_cnv_event.txt | awk '{print $2}'`; do echo ${i} > ${i}_interval.txt;done


                # Print the number of GAINS and LOSSES that DRAGEN found
                cat temp.vcf | grep -v '#' | awk '$7 == "PASS" {print}' | awk '{print $3}' | awk 'BEGIN {FS=":"} {print $2}' | sort | uniq -c | sort
            >>>
            output {
                File uncompressed_cnv_vcf = "temp.vcf"
                File cnv_events_file = "pass_cnv_event.txt"
                Array[File] cnv_intervals_files = glob("*_interval.txt")

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
    task length_check {
        input {
            File cnv_interval
        }
        command <<<
            set -e
            interval=$(cat ~{cnv_interval})

            start_pos=$(cat ~{cnv_interval} | cut -d ':' -f 2 | cut -d '-' -f 1)
            end_pos=$(cat ~{cnv_interval} | cut -d ':' -f 2 | cut -d '-' -f 2)

            event_length=$((end_pos - start_pos))
            echo $event_length > cnv_length.txt
            >>>
        output {
            Int cnv_event_length = read_int('cnv_length.txt')
        }
        runtime {
            docker: "us.gcr.io/broad-dsde-methods/liquidbiopsy:0.0.3.7"
            bootDiskSizeGb: 12
            cpu: 1
            memory: "8 GB"
            disks: "local-disk 100 HDD"
            preemptible: 2
            maxRetries: 3
        }
    }

    task blastn {
        input {
                String docker
                Int cnv_length
                File cnv_interval
                File cnv_vcf
                File reference_fasta
                File reference_blast_database
                File T2T_blast_database
            }
        command <<<
                set -e

                # Extract blast database from tar files
                mkdir -p /blastdb/reference_database
                mkdir -p /blastdb/t2t_database
                tar -xvf ~{reference_blast_database} -C /blastdb/reference_database/
                tar -xvf ~{T2T_blast_database} -C /blastdb/t2t_database/

                # Basename for the blast database
                reference_db_path=$(echo "`readlink -f /blastdb/reference_database/*`/`basename /blastdb/reference_database/*/*.nhr | cut -d '.' -f 1`")
                t2t_db_path=$(echo "`readlink -f /blastdb/t2t_database/*`/`basename /blastdb/t2t_database/*/*.nhr | cut -d '.' -f 1`")

                # Collect interval name
                interval=$(cat ~{cnv_interval})

                if ((~{cnv_length} > 15000)); then
                    cpu_num=16
                else
                    cpu_num=4
                fi

                # Try running Blastn
                mkdir /blastn/output
                (
                    python3 /blastn/main.py -i $interval -r ~{reference_fasta} -rd $reference_db_path -td $t2t_db_path -o /blastn/output -c $cpu_num
                ) || (
                    echo "$interval\tFail\tFail" > /blastn/output/${interval}_copy_number.txt
                )
                # Extract the CNV event name from the input file name
                filename=$(basename ~{cnv_vcf})
                # Unzip the input file if it is gzipped
                if test "${filename##*.}" = "gz"; then
                  gunzip -c ~{cnv_vcf} > temp.vcf
                else
                  cp ~{cnv_vcf} temp.vcf
                fi

                cat temp.vcf | grep $interval | awk '{print $5}' > cnv_event_type.txt
                paste cnv_event_type.txt /blastn/output/${interval}_copy_number.txt > annoed_blast_interval.txt
            >>>
        output {
            File blasted_cnv_interval = "annoed_blast_interval.txt"
        }
        runtime {
                docker: docker
                bootDiskSizeGb: 12
                cpu: if cnv_length > 15000 then 16 else 4
                memory: if cnv_length>15000 then "128GB" else "32GB"
                disks: "local-disk 100 HDD"
                preemptible: 2
                maxRetries: 3
        }
    }

    task gather_results {
        input{
            Array[File] cnv_blasted_results
         }
        command <<<
        set -e

        cat ~{sep=" " cnv_blasted_results} >> cnv_blasted.txt
            >>>
        output {
            File gathered_intervals_txt = "cnv_blasted.txt"
        }
        runtime {
                docker: "us.gcr.io/broad-dsde-methods/liquidbiopsy:0.0.3.7"
                bootDiskSizeGb: 12
                cpu: 1
                memory: "4 GB"
                disks: "local-disk 100 HDD"
                preemptible: 2
                maxRetries: 3
        }
    }




