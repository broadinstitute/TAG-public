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
                Array[File] cnv_intervals_chucks = glob("pass_cnv_chunk_*")

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