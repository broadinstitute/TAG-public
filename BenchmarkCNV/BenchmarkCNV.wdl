version 1.0

workflow Benchmark_CNV_Caller {
    input {
        String truth_sample_name
        String query_sample_name
        File truth_vcf
        File eval_cnv_vcf
        File eval_sv_vcf
        File wittyer_cnv_config
        File wittyer_sv_config
        File? bedfile
        String wittyer_cnv_evaluation_mode
        String wittyer_sv_evaluation_mode
        String wittyer_docker
        String wittyer4mat_docker
    }

    # benchmark cnv.vcf and sv.vcf using witty.er tool
    call BenchmarkCNV {
        input:
            wittyer_docker = wittyer_docker,
            truth_vcf = truth_vcf,
            truth_sample_name = truth_sample_name,
            query_sample_name = query_sample_name,
            eval_cnv_vcf = eval_cnv_vcf,
            cnv_config_file = wittyer_cnv_config,
            cnv_evaluation_mode = wittyer_cnv_evaluation_mode,
            eval_sv_vcf = eval_sv_vcf,
            sv_config_file = wittyer_sv_config,
            sv_evaluation_mode = wittyer_sv_evaluation_mode,
            bedfile = bedfile
    }

    # wittyer4mat to parse the wittyer json output
    call Wittyer4Mat {
        input:
            wittyer4mat_docker = wittyer4mat_docker,
            cnv_wittyer_stats = BenchmarkCNV.cnv_wittyer_stats,
            sv_wittyer_stats = BenchmarkCNV.sv_wittyer_stats,
            truth_sample_name = truth_sample_name
    }

    # Outputs that will be retained when execution is complete
    output {
        File cnv_wittyer_stats = BenchmarkCNV.cnv_wittyer_stats
        File cnv_wittyer_annotated_vcf = BenchmarkCNV.cnv_wittyer_annotated_vcf
        File cnv_wittyer_annotated_vcf_index = BenchmarkCNV.cnv_wittyer_annotated_vcf_index
        File sv_wittyer_stats = BenchmarkCNV.sv_wittyer_stats
        File sv_wittyer_annotated_vcf = BenchmarkCNV.sv_wittyer_annotated_vcf
        File sv_wittyer_annotated_vcf_index = BenchmarkCNV.sv_wittyer_annotated_vcf_index
        File cnv_deletion_stat = Wittyer4Mat.cnv_deletion_stat
        File cnv_duplication_stat = Wittyer4Mat.cnv_duplication_stat
        File sv_deletion_stat = Wittyer4Mat.sv_deletion_stat
        File sv_duplication_stat = Wittyer4Mat.sv_duplication_stat
        File sv_insertion_stat = Wittyer4Mat.sv_insertion_stat
    }
    meta {
        author: "Yueyao Gao"
        email: "gaoyueya@broadinstitute.org"
        description: "BenchmarkCNV.wdl is designed to evaluate the performance of Dragen CNV (Copy Number Variation) caller against GATK SV (Structural Variation) caller."
    }
}


    # Task 1: Benchmark the large variant vcf against truth set
    # If you are extracting vcf from a large callset vcf
    # Checkout /BenchmarkCNV/SelectSampleFromCallSet.wdl
    task BenchmarkCNV {

        input {
            String wittyer_docker
            File truth_vcf
            File eval_cnv_vcf
            File cnv_config_file
            String cnv_evaluation_mode
            File eval_sv_vcf
            File sv_config_file
            File? bedfile
            String sv_evaluation_mode
            String truth_sample_name
            String query_sample_name
            Int? mem
            Int? disk_space
            # If mem and disk size were not specified, use 4GB and 100 GB as default
            Int mem_size = select_first([mem, 4])
            Int disk_size = select_first([disk_space,100])
    }
        command <<<
            set -e

            if [[ -f "~{bedfile}" ]]; then
            # Run Benchmarking tool wittyer on dragen generated cnv.vcf with bed file
            /opt/Wittyer/Wittyer -i ~{eval_cnv_vcf} \
            -t ~{truth_vcf} \
            -em ~{cnv_evaluation_mode} \
            --configFile ~{cnv_config_file} \
            --includeBed ~{bedfile} \
            -o ~{truth_sample_name}_cnv_wittyer_output

            # Run Benchmarking tool wittyer on dragen generated sv.vcf with bed file
            /opt/Wittyer/Wittyer -i ~{eval_sv_vcf} \
            -t ~{truth_vcf} \
            -em ~{sv_evaluation_mode} \
            --configFile ~{sv_config_file} \
            --includeBed ~{bedfile} \
            -o ~{truth_sample_name}_sv_wittyer_output
            else
            # Run Benchmarking tool wittyer on dragen generated cnv.vcf
            /opt/Wittyer/Wittyer -i ~{eval_cnv_vcf} \
            -t ~{truth_vcf} \
            -em ~{cnv_evaluation_mode} \
            --configFile ~{cnv_config_file} \
            -o ~{truth_sample_name}_cnv_wittyer_output

            # Run Benchmarking tool wittyer on dragen generated sv.vcf
            /opt/Wittyer/Wittyer -i ~{eval_sv_vcf} \
            -t ~{truth_vcf} \
            -em ~{sv_evaluation_mode} \
            --configFile ~{sv_config_file} \
            -o ~{truth_sample_name}_sv_wittyer_output

        >>>
        runtime {
            docker: wittyer_docker
            bootDiskSizeGb: 12
            memory: mem_size + " GB"
            disks: "local-disk " + disk_size + " HDD"
            preemptible: 2
        }
        output {
            File cnv_wittyer_stats = "~{truth_sample_name}_cnv_wittyer_output/Wittyer.Stats.json"
            File cnv_wittyer_annotated_vcf = "~{truth_sample_name}_cnv_wittyer_output/Wittyer.~{truth_sample_name}.Vs.~{query_sample_name}.vcf.gz"
            File cnv_wittyer_annotated_vcf_index = "~{truth_sample_name}_cnv_wittyer_output/Wittyer.~{truth_sample_name}.Vs.~{query_sample_name}.vcf.gz.tbi"
            File sv_wittyer_stats = "~{truth_sample_name}_sv_wittyer_output/Wittyer.Stats.json"
            File sv_wittyer_annotated_vcf = "~{truth_sample_name}_sv_wittyer_output/Wittyer.~{truth_sample_name}.Vs.~{query_sample_name}.vcf.gz"
            File sv_wittyer_annotated_vcf_index = "~{truth_sample_name}_sv_wittyer_output/Wittyer.~{truth_sample_name}.Vs.~{query_sample_name}.vcf.gz.tbi"
        }
}
    # Task3: Format the wittyer json output
    task Wittyer4Mat{

        input{
            String wittyer4mat_docker
            File cnv_wittyer_stats
            File sv_wittyer_stats
            String truth_sample_name
        }
        command <<<
            set -e

            # Run wittyer4mat script on cnv wittyer output
            mkdir ~{truth_sample_name}_cnv_wittyer4mat
            conda run --no-capture-output \
            -n wittyer-parser \
            python3 /wittyer4mat/wittyer_4mat.py -i ~{cnv_wittyer_stats} \
            -t cnv \
            -o ~{truth_sample_name}_cnv_wittyer4mat

            # Run wittyer4mat script on sv wittyer output
            mkdir ~{truth_sample_name}_sv_wittyer4mat
            conda run --no-capture-output \
            -n wittyer-parser \
            python3 /wittyer4mat/wittyer_4mat.py -i ~{sv_wittyer_stats} \
            -t sv \
            -o ~{truth_sample_name}_sv_wittyer4mat
    >>>
        runtime {
            docker: wittyer4mat_docker
            preemptible: 2
        }
        output {
            File cnv_deletion_stat = "~{truth_sample_name}_cnv_wittyer4mat/wittyer_cnv_Deletion_output.csv"
            File cnv_duplication_stat = "~{truth_sample_name}_cnv_wittyer4mat/wittyer_cnv_Duplication_output.csv"
            File sv_deletion_stat = "~{truth_sample_name}_sv_wittyer4mat/wittyer_sv_Deletion_output.csv"
            File sv_duplication_stat = "~{truth_sample_name}_sv_wittyer4mat/wittyer_sv_Duplication_output.csv"
            File sv_insertion_stat = "~{truth_sample_name}_sv_wittyer4mat/wittyer_sv_Insertion_output.csv"
        }
    }