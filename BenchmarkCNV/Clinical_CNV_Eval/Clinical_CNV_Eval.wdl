version 1.0

workflow Clinical_CNV_Eval {
    input {
        String Truth_Type
        String Truth_Coordinates
        File Truth_VCF_Header
        File Pre_Wittyer_Config
        Int Wittyer_BPD
        Float Wittyer_PD
        File Sample_VCF
        String Wittyer_Docker = "yg96/wittyer:v2"
        File Proband_ShortVariant_VCF
        File Normal_ShortVariant_VCF
    }
    call GenerateTruthVCF {
        input:
            Truth_Type = Truth_Type,
            Truth_Coordinates = Truth_Coordinates,
            Truth_VCF_Header = Truth_VCF_Header
    }
    call SetWittyerConfig {
        input:
            Pre_Wittyer_Config = Pre_Wittyer_Config,
            Wittyer_BPD = 750000,
            Wittyer_PD = 0.25
    }
    call CompareVCF {
        input:
            Truth_VCF = GenerateTruthVCF.truth_vcf,
            Sample_VCF = Sample_VCF,
            wittyer_config = SetWittyerConfig.wittyer_config,
            Wittyer_Docker = Wittyer_Docker
    }
    call PostProcessWittyer {
        input:
            Comparison_VCF = CompareVCF.comparison_vcf,
            Wittyer_BPD = Wittyer_BPD,
            Wittyer_PD = Wittyer_PD,
            Proband_ShortVariant_VCF = Proband_ShortVariant_VCF,
            Normal_ShortVariant_VCF = Normal_ShortVariant_VCF
    }
    output {
        File Truth_VCF = GenerateTruthVCF.truth_vcf
        File Comparison_VCF = CompareVCF.comparison_vcf
        Int truth_length = PostProcessWittyer.truth_length
        String wittyer_decision = PostProcessWittyer.wittyer_decision
        File overlapping_dragen_calls = PostProcessWittyer.overlapping_dragen_calls
        File constructed_dragen_call = PostProcessWittyer.constructed_dragen_call
        String merged_eval_decision = PostProcessWittyer.merged_eval_decision
        String merged_overlap_ratio = PostProcessWittyer.merged_overlap_ratio
        File? comparison_plot = PostProcessWittyer.comparison_plot
        File proband_snp_af_within_truth_interval_stats = PostProcessWittyer.proband_snp_af_within_truth_interval_stats
        File normal_snp_af_within_truth_interval_stats = PostProcessWittyer.normal_snp_af_within_truth_interval_stats
        File heterozyogisty_within_truth_interval_plot = PostProcessWittyer.heterozyogisty_within_truth_interval_plot
    }
}
    task GenerateTruthVCF{
        input {
            String Truth_Type
            String Truth_Coordinates
            File Truth_VCF_Header
        }
        command <<<
python <<CODE

cnvtype_dict = {"Gain":"DUP","Loss":"DEL"}

with open('truth.vcf','w') as vcf_file, open("~{Truth_VCF_Header}",'r') as vcf_header:
    # write header
    for line in vcf_header:
        vcf_file.write(line)
    # write truth variant
    chrom = "~{Truth_Coordinates}".split(':')[0].split('chr')[1]
    pos = "~{Truth_Coordinates}".split(':')[1].split('-')[0]
    end = "~{Truth_Coordinates}".split(':')[1].split('-')[1]
    svlen = abs(int(pos)-int(end))+1
    cnv_type = cnvtype_dict["~{Truth_Type}"]

    vcf_file.write(f"{chrom}\t{pos}\t~{Truth_Coordinates}\tN\t<{cnv_type}>\t.\tPASS\tEND={end};SVTYPE={cnv_type};SVLEN={svlen}\t.\t.\n")
CODE
    >>>
        runtime {
            docker: "us.gcr.io/tag-team-160914/wittyer4mat:v13-pure-env"
            preemptible: 2
        }
        output {
            File truth_vcf = "truth.vcf"
        }
    }
    task SetWittyerConfig{
        input {
            File Pre_Wittyer_Config
            Int Wittyer_BPD
            Float Wittyer_PD
        }
        command <<<
    set -e
    # Set wittyer config with input parameters
    cat ~{Pre_Wittyer_Config} | jq '[.[] | .bpDistance = ~{Wittyer_BPD} | .percentDistance = ~{Wittyer_PD}]' > wittyer_config.json

    >>>
        runtime {
            docker: "stedolan/jq:latest"
            preemptible: 2
        }
        output {
            File wittyer_config = "wittyer_config.json"
        }
}
    task CompareVCF{
        input {
            File Truth_VCF
            File Sample_VCF
            File wittyer_config
            String Wittyer_Docker
        }
        command <<<
    set -e
    # Run Benchmarking tool wittyer on dragen generated cnv_sv.vcf and clinical truth vcf
    /opt/Wittyer/Wittyer -i ~{Sample_VCF} \
    -t ~{Truth_VCF} \
    -em CrossTypeAndSimpleCounting \
    --configFile ~{wittyer_config} \
    -o wittyer_output

    # Move wittyer output to the output directory
    mv wittyer_output/Wittyer.*.Vs.*.vcf.gz comparison.vcf.gz
    mv wittyer_output/Wittyer.*.Vs.*.vcf.gz.tbi comparison.vcf.gz.tbi
    >>>
        runtime {
            docker: Wittyer_Docker
            preemptible: 2
        }
        output {
            File comparison_vcf = "comparison.vcf.gz"
            File comparison_vcf_tbi = "comparison.vcf.gz.tbi"
            File wittyer_stats = "wittyer_output/Wittyer.Stats.json"
        }
    }
    task PostProcessWittyer{
        input {
            File Comparison_VCF
            Int Wittyer_BPD
            Float Wittyer_PD
            String Wittyer_Postprocess_Docker = "us.gcr.io/tag-team-160914/cnv-clinical-eval:0.1.1-alpha"
            File Proband_ShortVariant_VCF
            File Normal_ShortVariant_VCF
            Int memory = 16
            Int cpu = 4
        }
        command <<<
    set -e
            # post process wittyer output
            # Check TP, FN and see if merge is needed
            conda run --no-capture-output \
            -n Clinical-CNV-Env \
            python3 /BaseImage/CNV_Clinical_Eval/scripts/cnv_clinical_eval.py \
            --vcf ~{Comparison_VCF} \
            --bpd ~{Wittyer_BPD} \
            --pd ~{Wittyer_PD} \
            --proband_vcf ~{Proband_ShortVariant_VCF} \
            --normal_vcf ~{Normal_ShortVariant_VCF}

            # rename the DRAGEN-vs-Truth output figure if it exists
            if [ -f *_compare_sv.png ]; then
                mv *_compare_sv.png wittyer_comparison_plot.png
            fi

            # Rename the truth ROH status file
            mv *_AF_plots.png proband_vs_normal_af_plot.png



    >>>
        runtime {
            docker: Wittyer_Postprocess_Docker
            preemptible: 2
            memory: memory + " GB"
            cpu: cpu
            disks: "local-disk 100 HDD"
        }
        output {
            Int truth_length = read_int("truth_length.txt")
            String wittyer_decision = read_string("wittyer_decision.txt")
            File overlapping_dragen_calls = "dragen_calls.txt"
            File constructed_dragen_call = "constructed_interval.txt"
            String merged_eval_decision = read_string("merged_eval_decision.txt")
            String merged_overlap_ratio = read_string("merged_overlap_ratio.txt")
            File? comparison_plot = "wittyer_comparison_plot.png"
            File proband_snp_af_within_truth_interval_stats = "proband_snp_within_interval_af.csv"
            File normal_snp_af_within_truth_interval_stats = "normal_snp_within_interval_af.csv"
            File heterozyogisty_within_truth_interval_plot = "proband_vs_normal_af_plot.png"
        }
    }
    