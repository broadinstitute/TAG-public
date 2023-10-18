version 1.0

workflow Benchmark_CNV_Caller {
    input {
        String truth_sample_name
        String query_sample_name
        File truth_vcf
        File eval_vcf
        File wittyer_config
        File? bedfile
        String wittyer_evaluation_mode
        String wittyer_docker = "yg96/wittyer:v2"
        String wittyer4mat_docker = "us.gcr.io/tag-team-160914/wittyer4mat:v12-prestat-broom"
        Boolean run_wittyer4mat
        Boolean get_queryfp
    }

    # benchmark cnv.vcf and sv.vcf using witty.er tool
    call BenchmarkCNV {
        input:
            wittyer_docker = wittyer_docker,
            truth_vcf = truth_vcf,
            truth_sample_name = truth_sample_name,
            query_sample_name = query_sample_name,
            eval_vcf = eval_vcf,
            wittyer_config = wittyer_config,
            wittyer_evaluation_mode = wittyer_evaluation_mode,
            bedfile = bedfile
    }

    # wittyer4mat to parse the wittyer json output
    if (run_wittyer4mat) {
        call Wittyer4Mat {
            input:
                wittyer4mat_docker = wittyer4mat_docker,
                wittyer_stats = BenchmarkCNV.wittyer_stats,
                truth_sample_name = truth_sample_name
        }
    }

    # Grab Query FP to the datatable
        if (get_queryfp) {
        call queryFP {
            input:
                wittyer4mat_docker = wittyer4mat_docker,
                wittyer_stats = BenchmarkCNV.wittyer_stats
        }
    }

    # Outputs that will be retained when execution is complete
    output {
        File wittyer_stats = BenchmarkCNV.wittyer_stats
        File wittyer_annotated_vcf = BenchmarkCNV.wittyer_annotated_vcf
        File wittyer_annotated_vcf_index = BenchmarkCNV.wittyer_annotated_vcf_index
        Array[File]? Wittyer4Mat_event_stats = Wittyer4Mat.event_level_wittyer_stats
        Int? queryFP = queryFP.queryFP
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
            File eval_vcf
            File wittyer_config
            String wittyer_evaluation_mode
            File? bedfile
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

            # Run Benchmarking tool wittyer on dragen generated cnv.vcf and truth set
            /opt/Wittyer/Wittyer -i ~{eval_vcf} \
            -t ~{truth_vcf} \
            -em ~{wittyer_evaluation_mode} \
            --configFile ~{wittyer_config} \
            ~{'--includeBed '+ bedfile} \
            -o ~{truth_sample_name}_wittyer_output

            # Move wittyer output to the output directory
            mv ~{truth_sample_name}_wittyer_output/Wittyer.*.Vs.*.vcf.gz ~{truth_sample_name}_wittyer_output/Wittyer.~{truth_sample_name}.Vs.~{query_sample_name}.vcf.gz
            mv ~{truth_sample_name}_wittyer_output/Wittyer.*.Vs.*.vcf.gz.tbi ~{truth_sample_name}_wittyer_output/Wittyer.~{truth_sample_name}.Vs.~{query_sample_name}.vcf.gz.tbi

        >>>
        runtime {
            docker: wittyer_docker
            bootDiskSizeGb: 12
            memory: mem_size + " GB"
            disks: "local-disk " + disk_size + " HDD"
            preemptible: 2
        }
        output {
            File wittyer_stats = "~{truth_sample_name}_wittyer_output/Wittyer.Stats.json"
            File wittyer_annotated_vcf = "~{truth_sample_name}_wittyer_output/Wittyer.~{truth_sample_name}.Vs.~{query_sample_name}.vcf.gz"
            File wittyer_annotated_vcf_index = "~{truth_sample_name}_wittyer_output/Wittyer.~{truth_sample_name}.Vs.~{query_sample_name}.vcf.gz.tbi"
        }
}
    # Task2: Format the wittyer json output
    task Wittyer4Mat{

        input{
            String wittyer4mat_docker
            File wittyer_stats
            String truth_sample_name
        }
        command <<<
            set -e

            # Run wittyer4mat script on wittyer output
            conda run --no-capture-output \
            -n wittyer-parser \
            python3 /wittyer4mat/wittyer_4mat.py -i ~{wittyer_stats} \
            -t event \
            -o ~{truth_sample_name}_event_level_wittyer4mat
    >>>
        runtime {
            docker: wittyer4mat_docker
            preemptible: 2
        }
        output {
            Array[File] event_level_wittyer_stats = glob("~{truth_sample_name}_event_level_wittyer4mat/*.csv")
        }
    }

    #Task3: Grab Query FP to the datatable
    task queryFP{
        input{
            String wittyer4mat_docker
            File wittyer_stats
        }
        command <<<
        set -e
        pip install dictor
        python <<CODE
            import json
            import dictor
            with open("~{wittyer_stats}") as data:
                data = json.load(data)

            queryFP = dictor(data, f"PerSampleStats.0.OverallStats.0.QueryFpCount")
            with open("queryFP.txt", "w") as f:
                f.write(str(queryFP))
      CODE
    >>>
        runtime {
            docker: wittyer4mat_docker
            preemptible: 2
        }
        output {
            Int queryFP = read_int("queryFP.txt")
        }
    }
