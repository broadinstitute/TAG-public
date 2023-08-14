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
#                    reference_last_database = reference_last_database,
#                    T2T_last_database = T2T_last_database,
                    last_docker = last_docker
            }
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
#            File reference_last_database
#            File T2T_last_database
            String last_docker
        }
        command <<<
            set -e

#            # Extract blast database from tar files
#            mkdir -p /lastdb/reference_database
#            mkdir -p /lastdb/t2t_database
#            tar -xvf {reference_last_database} -C /lastdb/reference_database/
#            tar -xvf {T2T_last_database} -C /lastdb/t2t_database/
#            # Basename for the blast database
#            # If there are dot in the name, it will cause error in extracting the basename
#            reference_db_path=$(echo "`readlink -f /lastdb/reference_database/*`/`basename /lastdb/reference_database/*/*.prj | cut -d '.' -f 1`")
#            t2t_db_path=$(echo "`readlink -f /lastdb/t2t_database/*`/`basename /lastdb/t2t_database/*/*.prj | cut -d '.' -f 1`")

            for i in `ls ~{sep=" " interval_sequence_fasta}`; do
                interval_name=$(basename ${i} ".fasta")
                echo $interval_name
#                /last-1460/bin/lastal5 $reference_db_path {interval_sequence_fasta} -v -P 0 -l 30 -f BlastTab > ${interval_name}_ref_lastal_alignment.txt
#                /last-1460/bin/lastal5 $t2t_db_path {interval_sequence_fasta} -v -P 0 -l 30 -f BlastTab > ${interval_name}_t2t_lastal_alignment.txt
            done

    >>>

        runtime {
                docker: last_docker
                bootDiskSizeGb: 12
                cpu: 16
                memory: "128 GB"
                disks: "local-disk 500 HDD"
                preemptible: 0
                maxRetries: 0
        }
    }