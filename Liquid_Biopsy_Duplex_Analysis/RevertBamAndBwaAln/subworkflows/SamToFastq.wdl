workflow samToFastqTest {
	call samtofastq
}

task samtofastq {

    File input_bam

    Float? memory = 15
    Int? disk_space = 200
    Int? num_threads = 4
    Int num_preempt

    String? docker_override
    String gatk_docker = select_first([docker_override,"us.gcr.io/broad-gatk/gatk:4.5.0.0"])
    String? gatk_override

    command {
        set -euo pipefail
        
        export GATK_PATH=${default="/gatk/gatk" gatk_override}

        mkdir -p samtofastq  # workaround for named pipes
        
        samtools view -H ${input_bam} | grep ^@RG > read_groups.txt

        $GATK_PATH SamToFastq \
        -I ${input_bam} \
        --OUTPUT_PER_RG true \
        --COMPRESS_OUTPUTS_PER_RG true \
        --VALIDATION_STRINGENCY SILENT \
        --OUTPUT_DIR samtofastq 
        
        mv samtofastq/*.fastq.gz .
    }

    output {
        Array[File] firstEndFastqs = glob("*1.fastq.gz")
        Array[File] secondEndFastqs = glob("*2.fastq.gz")
        Array[String] read_groups = read_lines("read_groups.txt")
    }

    runtime {
        docker: gatk_docker
        memory: select_first([memory])+"GB"
        disks: "local-disk "+ select_first([disk_space])+" HDD"
        cpu: select_first([num_threads])
        preemptible: "${num_preempt}"
    }
}