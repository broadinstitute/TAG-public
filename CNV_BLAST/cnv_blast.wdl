version 1.0

    workflow cnv_blast{
        input {
            File cnv_vcf
        }
        call extract_cnv {
            input:
                cnv_vcf = cnv_vcf
        }

        scatter (cnv_event_file in extract_cnv.cnv_intervals_files) {

            call test_task {
                input:
                    cnv_interval = cnv_event_file
            }
        }

        call gather_results {
            input:
                cnv_interval_txt = test_task.cnv_interval_txt
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
    task test_task {
        input {
            File cnv_interval
            }
        command <<<
            set -e

            echo `cat ~{cnv_interval}` + " is a CNV interval" > ~{cnv_interval}.txt
            >>>
        output {
            File cnv_interval_txt = "~{cnv_interval}.txt"
        }
        runtime {
                docker: "us.gcr.io/broad-dsde-methods/liquidbiopsy:0.0.3.7"
                bootDiskSizeGb: 12
                cpu: 1
                memory: "4 GB"
                disks: "local-disk 10 HDD"
                preemptible: 2
                #maxRetries: 3
        }
    }

    task gather_results {
        input{
            Array[File] cnv_interval_txt
         }
        command <<<
        set -e

        cat ~{sep=" " cnv_interval_txt} >> cnv_intervals.txt
            >>>
        output {
            File gathered_intervals_txt = "cnv_intervals.txt"
        }
        runtime {
                docker: "us.gcr.io/broad-dsde-methods/liquidbiopsy:0.0.3.7"
                bootDiskSizeGb: 12
                cpu: 1
                memory: "4 GB"
                disks: "local-disk 10 HDD"
                preemptible: 2
                #maxRetries: 3
        }
    }




