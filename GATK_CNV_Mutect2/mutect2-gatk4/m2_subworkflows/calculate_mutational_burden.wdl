version 1.0

workflow calculate_mutational_burden {
    input{
    File ref_fasta
    File ref_fai
    File ref_dict
    File tumor_reads
    File tumor_reads_index
    File? normal_reads
    File? normal_reads_index

    File input_maf
    File? intervals
    File mutburden_python_script
    File context_script_override
    String tag_docker
    String output_basename

    Int disk_pad
    Int tumor_size
    Int normal_size
    Int ref_size

    Int? preemptible_attempts
    }

    if (basename(tumor_reads) != basename(tumor_reads, ".cram")) {
        call CramToBam as TumorCramToBam {
            input:
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                cram = tumor_reads,
                crai = tumor_reads_index,
                name = output_basename,
                disk_size = tumor_size + normal_size + ref_size + disk_pad
        }
    }

    String normal_or_empty = select_first([normal_reads, ""])
    if (basename(normal_or_empty) != basename(normal_or_empty, ".cram")) {
        String normal_basename = basename(basename(normal_or_empty, ".bam"),".cram")
        call CramToBam as NormalCramToBam {
            input:
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                cram = normal_reads,
                crai = normal_reads_index,
                name = normal_basename,
                disk_size = tumor_size + normal_size + ref_size + disk_pad
        }
    }

    File tumor_bam = select_first([TumorCramToBam.output_bam, tumor_reads])
    File tumor_bai = select_first([TumorCramToBam.output_bai, tumor_reads_index])
    File? normal_bam = if defined(normal_reads) then select_first([NormalCramToBam.output_bam, normal_reads]) else normal_reads
    File? normal_bai = if defined(normal_reads) then select_first([NormalCramToBam.output_bai, normal_reads_index]) else normal_reads_index

    Int tumor_bam_size = ceil(size(tumor_bam, "GB") + size(tumor_bai, "GB"))
    Int normal_bam_size = if defined(normal_bam) then ceil(size(normal_bam, "GB") + size(normal_bai, "GB")) else 0

    call CallableLoci {
        input:
            output_basename = output_basename,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            tumor_bam = tumor_bam,
            tumor_bai = tumor_bai,
            normal_bam = normal_bam,
            normal_bai = normal_bai,
            intervals = intervals,
            tag_docker = tag_docker,
            context_script_override = context_script_override,
            preemptible_attempts = preemptible_attempts,
            disk_space = tumor_size + normal_size + ref_size + disk_pad
    }

    call MutationalBurden {
        input:
            output_basename = output_basename,
            input_maf = input_maf,
            mutburden_python_script=mutburden_python_script,
            tag_docker = tag_docker,
            callable_bases = CallableLoci.callable_bases,
            preemptible_attempts = preemptible_attempts,
            disk_space = ceil(size(input_maf, "GB") + disk_pad)
    }

    output {
        String callable_bases = CallableLoci.callable_bases
        File callable_regions = CallableLoci.callable_regions
        File callable_contexts = CallableLoci.callable_contexts

        String total_variants = MutationalBurden.total_variants
        String coding_variants = MutationalBurden.coding_variants
        String coding_mutations_per_mb = MutationalBurden.coding_mutations_per_mb
        String noncoding_variants = MutationalBurden.noncoding_variants
        String noncoding_mutations_per_mb = MutationalBurden.noncoding_mutations_per_mb
        File mutational_burden = MutationalBurden.mutational_burden
    }
}

task CallableLoci {
    input{
    String output_basename
    File ref_fasta
    File ref_fai
    File ref_dict
    File tumor_bam
    File tumor_bai
    File? normal_bam
    File? normal_bai
    File? intervals

    String tag_docker
    File? gatk_override
    File? context_script_override

    Int? preemptible_attempts
    Int? disk_space
    Int? mem
    Int? cpu

    # Cutoff to judge covered bases
    Int? tumor_coverage
    Int? normal_coverage

}

    Int machine_mem = select_first([mem, 4])*1000
    Int command_mem = machine_mem - 500
    Int tumor_cutoff = select_first([tumor_coverage,14])
    Int normal_cutoff = select_first([normal_coverage,8])

    command <<<
        set -e
        export GATK_JAR=~{default="/usr/tag/GATK36.jar" gatk_override}
        export CONTEXT_PY=~{default="/usr/tag/kmer_freq.py" context_script_override}

        java "-Xmx~{command_mem}m" -jar $GATK_JAR -T CallableLoci \
            -I ~{tumor_bam} \
            -R ~{ref_fasta} \
            --minMappingQuality 20 \
            --minBaseQuality 20 \
            --minDepth ~{tumor_cutoff} \
            ~{"-L " + intervals} \
            -o tumor_callable.bed \
            --summary tumor_callable.summary

        if [[ -f "~{normal_bam}" ]]; then
            java "-Xmx~{command_mem}m" -jar $GATK_JAR -T CallableLoci \
                -I ~{normal_bam} \
                -R ~{ref_fasta} \
                --minMappingQuality 20 \
                --minBaseQuality 20 \
                --minDepth ~{normal_cutoff} \
                ~{"-L " + intervals} \
                -o normal_callable.bed \
                --summary normal_callable.summary

            bedtools intersect -a <(grep 'CALLABLE' tumor_callable.bed) \
                               -b <(grep 'CALLABLE' normal_callable.bed) > ~{output_basename}_callable.bed
        else
            grep 'CALLABLE' tumor_callable.bed > ~{output_basename}_callable.bed
        fi

        # Tally callable bases from BED
        if [[ -s "~{output_basename}_callable.bed" ]]; then
            awk 'BEGIN{sum=0}{sum+=$3-$2}END{print(sum)}' ~{output_basename}_callable.bed > callable_bases.txt
        else
            echo "0" > callable_bases.txt
        fi

        # Obtain callable bases in 3-base contexts
        # awk command is for including flanking bases
        awk 'BEGIN{OFS="\t"; FS="\t"}{$2-=1; $3+=1; print $0}' ~{output_basename}_callable.bed | \
        bedtools getfasta -fi ~{ref_fasta} -bed stdin | \
        python ${CONTEXT_PY} 3 - > ~{output_basename}_context.txt
    >>>

    runtime {
        docker: tag_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 12]) + " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
    }

    output {
        String callable_bases = read_string("callable_bases.txt")
        File callable_regions = "${output_basename}_callable.bed"
        File callable_contexts = "${output_basename}_context.txt"
    }
}

task MutationalBurden {
    input{
    String output_basename
    File input_maf
    String callable_bases
    File mutburden_python_script

    # runtime
    String tag_docker
    Int? preemptible_attempts
    Int? disk_space
    Int? mem


}

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    command <<<

        python ~{mutburden_python_script} --sample-id ~{output_basename} ~{callable_bases} ~{input_maf}

        # Extract values for displaying in the Terra/FireCloud data table
        grep "^total_variants" ~{output_basename}.mutational_burden.txt | cut -f2 > total_variants.txt
        grep "^coding_variants" ~{output_basename}.mutational_burden.txt | cut -f2 > coding_variants.txt
        grep "^noncoding_variants" ~{output_basename}.mutational_burden.txt | cut -f2 > noncoding_variants.txt
        grep "^coding_mutations_per_Mb" ~{output_basename}.mutational_burden.txt | cut -f2 > coding_mb.txt
        grep "^noncoding_mutations_per_Mb" ~{output_basename}.mutational_burden.txt | cut -f2 > noncoding_mb.txt
    >>>

    output  {
        File mutational_burden="~{output_basename}.mutational_burden.txt"
        String total_variants = read_string("total_variants.txt")
        String coding_variants = read_string("coding_variants.txt")
        String noncoding_variants = read_string("noncoding_variants.txt")
        String coding_mutations_per_mb = read_string("coding_mb.txt")
        String noncoding_mutations_per_mb = read_string("noncoding_mb.txt")
    }

    runtime {
        docker: tag_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 10]) + " HDD"
        preemptible: select_first([preemptible_attempts, 10])
    }
}

task CramToBam {
    input {
      File ref_fasta
      File ref_fai
      File ref_dict
      #cram and crai must be optional since Normal cram is optional
      File? cram
      File? crai
      String name
      Int disk_size
      Int? mem
    }

    Int machine_mem = if defined(mem) then mem * 1000 else 6000

    #Calls samtools view to do the conversion
    command {
        #Set -e and -o says if any command I run fails in this script, make sure to return a failure
        set -e
        set -o pipefail

        samtools view -h -T ~{ref_fasta} ~{cram} |
            samtools view -b -o ~{name}.bam -
        samtools index -b ~{name}.bam
        mv ~{name}.bam.bai ~{name}.bai
    }

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
        memory: machine_mem + " MB"
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File output_bam = "~{name}.bam"
        File output_bai = "~{name}.bai"
    }
}