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
            Wittyer_BPD = 250000,
            Wittyer_PD = 0.25
    }
    call CompareVCF {
        input:
            Truth_VCF = GenerateTruthVCF.truth_vcf,
            Sample_VCF = Sample_VCF,
            wittyer_config = SetWittyerConfig.wittyer_config,
            Wittyer_Docker = Wittyer_Docker
    }
    call CheckTPVariant {
        input:
            Comparison_VCF = CompareVCF.comparison_vcf
    }
    output {
        File Truth_VCF = GenerateTruthVCF.truth_vcf
        File Comparison_VCF = CompareVCF.comparison_vcf
        Int Called = CheckTPVariant.tp_called
        String Called_coordinates = CheckTPVariant.eval_coordinates
        Int Called_length = CheckTPVariant.eval_length
        String Called_svclaim = CheckTPVariant.eval_svclaim
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
            File comparison_vcf = "compare.vcf.gz"
            File comparison_vcf_tbi = "compare.vcf.gz.tbi"
            File wittyer_stats = "wittyer_output/Wittyer.Stats.json"
        }
    }
    task CheckTPVariant{
        input {
            File Comparison_VCF
        }
        command <<<
    set -e
    # Check if TP variant is called
    # If TP variant is called, return 1, else return 0
    # I hardcoded DRAGEN here. But it can be changed to any other caller
    if zcat ~{Comparison_VCF} | grep TP | grep DRAGEN | awk '{print $3}' | grep -q .; then
        echo "1" > called.txt
    else
        echo "0" > called.txt
    fi

    zcat ~{Comparison_VCF} | grep TP | grep DRAGEN | awk '{print $3}' > eval_coordinates.txt
    zcat ~{Comparison_VCF} | grep TP | grep DRAGEN  | awk '{print $8}' | grep -o 'SVLEN=[^;]*' | awk -F'=' '{print $2}' > eval_svlength.txt
    zcat ~{Comparison_VCF} | grep TP | grep DRAGEN  | awk '{print $8}' | grep -o 'SVCLAIM=[^;]*' | awk -F'=' '{print $2}' > eval_svclaim.txt
    >>>
        runtime {
            docker: "ubuntu:latest"
            preemptible: 2
        }
        output {
            Int tp_called = read_int("called.txt")
            String eval_coordinates = read_string("eval_coordinates.txt")
            Int eval_length = read_int("eval_svlength.txt")
            Int eval_svclaim = read_string("eval_svclaim.txt")
        }
    }
    