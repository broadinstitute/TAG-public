version 1.0

# Main workflow: performs evaluation of query_vcf against base_vcf using truvari bench
workflow BenchmarkSV{
    input {
        File base_vcf
        File base_vcf_index
        String base_vcf_sample_name

        File query_vcf
        File query_vcf_index
        String query_vcf_sample_name

        File evaluation_intervals
    }
    call truvari_bench {
        input:
            base_vcf = base_vcf,
            base_vcf_index = base_vcf_index,
            base_vcf_sample_name = base_vcf_sample_name,
            query_vcf = query_vcf,
            query_vcf_index = query_vcf_index,
            query_vcf_sample_name = query_vcf_sample_name,
            evaluation_intervals = evaluation_intervals
    }
    output {
        File bench_summary = truvari_bench.bench_summary
        File tp_base_vcf = truvari_bench.tp_base_vcf
        File tp_comp_vcf = truvari_bench.tp_comp_vcf
        File fp_vcf = truvari_bench.fp_vcf
        File fn_vcf = truvari_bench.fn_vcf
    }
    meta {
        author: "Yueyao Gao"
        email: "gaoyueya@broadinstitute.org"
        description: "BenchmarkSV.wdl is designed to run truvari bench on a single sample. The input files are a base vcf file, a query vcf file, and a bed file with evaluation intervals. The output files are truvari bench results and a sheet with benchmark statistics stratified by SVTYPE."
    }
}

task truvari_bench {
    input {
        File base_vcf
        File base_vcf_index
        String base_vcf_sample_name

        File query_vcf
        File query_vcf_index
        String query_vcf_sample_name

        File evaluation_intervals
        String truvari_docker = "us.gcr.io/tag-public/truvari:v5.0.0"
    }
    command <<<
        set -e
        mkdir truvari_output

        truvari bench \
        -b ~{base_vcf} \
        -c ~{query_vcf} \
        --dup-to-ins \
        --pctseq 0 \
        --refdist=2000 \
        --chunksize=2000 \
        --pctsize=0.7 \
        --pctovl=0.0 \
        --passonly \
        --minhaplen=50 \
        --sizemin=50 \
        --sizefilt=35 \
        --sizemax=500000000 \
        --no-ref=c \
        --pick=ac \
        --extend=0 \
        --includebed ~{evaluation_intervals} \
        -o truvari_output
    >>>
    output {
        File bench_summary = "truvari_output/summary.json"
        File tp_base_vcf = "truvari_output/tp-base.vcf.gz"
        File tp_base_vcf_index = "truvari_output/tp-base.vcf.gz.tbi"
        File tp_comp_vcf = "truvari_output/tp-comp.vcf.gz"
        File tp_comp_vcf_index = "truvari_output/tp-comp.vcf.gz.tbi"
        File fp_vcf = "truvari_output/fp.vcf.gz"
        File fp_vcf_index = "truvari_output/fp.vcf.gz.tbi"
        File fn_vcf = "truvari_output/fn.vcf.gz"
        File fn_vcf_index = "truvari_output/fn.vcf.gz.tbi"
    }
    runtime {
        docker: truvari_docker
        memory: "32 GB"
        cpu: 4
        disk_size_gb: 500
        preemptible: 3
    }
}