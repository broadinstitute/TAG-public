version 1.0

    workflow cnv_blast{
        input {
            File cnv_vcf
        }
        call extract_cnv {
            input:
                cnv_vcf = cnv_vcf
        }

        output {
            File cnv_events = extract_cnv.cnv_events
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
                # Extract the PASS events from the DRAGEN CNV VCF and save them to a file
                cat ~{cnv_vcf}| grep -v '#' | awk '$7 == "PASS" {print}' | awk '{print $3}' > pass_cnv_events.txt

                # Print the number of GAINS and LOSSES that DRAGEN found
                cat NA24385.cnv_sv.vcf| grep -v '#' | awk '$7 == "PASS" {print}' | awk '{print $3}' | awk 'BEGIN {FS=":"} {print $2}' | sort | uniq -c | sort
            >>>
            output {
                File cnv_events = "pass_cnv_events.txt"
            }
            runtime {
                docker: "us.gcr.io/broad-dsde-methods/liquidbiopsy:0.0.3.7"
                bootDiskSizeGb: 12
                cpu: cpu_num
                memory: mem_size + " GB"
                disks: "local-disk " + disk_size + " HDD"
                preemptible: 2
                maxRetries: 3
            }
        }


