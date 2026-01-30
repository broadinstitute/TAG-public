version 1.0

workflow FilterViralReads {
    input {
        File input_bam
        File viral_reference
        String sample_name
    }

    call FilterViralBam {
        input:
            bam_file = input_bam,
            reference_fasta = viral_reference,
            basename = sample_name
    }
    
    output {
        File viral_bam = FilterViralBam.out_bam
        File viral_bai = FilterViralBam.out_bai
    }
}

task FilterViralBam {
    input {
        File bam_file
        File reference_fasta
        String basename
        
        # Runtime attributes (can be overridden in inputs)
        Int threads = 8
        Int memory_gb = 16
        # Replace this string with the tag of the image you built in Part 1
        String docker_image = "fleharty/viral-bam-filter:v1" 
    }

    command <<<
        set -e
        
        echo "Processing Sample: ~{basename}"
        
        # 1. Convert BAM to FASTA
        # We use ~{threads} to utilize the requested CPU cores
        samtools fastq -@ ~{threads} "~{bam_file}" -o "~{basename}.fasta"

        # 2. Run BBDuk
        # Note: We set -Xmx to (memory_gb - 4) to leave room for overhead
        # We reference the WDL input ~{reference_fasta} directly
        bbduk.sh -Xmx~{memory_gb - 4}g \
            in="~{basename}.fasta" \
            ref="~{reference_fasta}" \
            outm="~{basename}.viral.fastq" \
            k=19 hdist=0 rcomp=t \
            threads=~{threads} tbo=f tpe=f

        # Cleanup intermediate FASTA
        rm "~{basename}.fasta"

        # 3. Extract read names
        # Standard grep/sed pipeline
        grep '^@' "~{basename}.viral.fastq" | sed 's/^@//' | cut -d ' ' -f1 | sort -u > "~{basename}.viral.names"

        # 4. Subset original BAM using the read names
        samtools view -N "~{basename}.viral.names" -b "~{bam_file}" > "~{basename}.viral.bam"
        samtools index "~{basename}.viral.bam"

        # Final Cleanup
        rm "~{basename}.viral.names"
        rm "~{basename}.viral.fastq"
    >>>

    output {
        File out_bam = "~{basename}.viral.bam"
        File out_bai = "~{basename}.viral.bam.bai"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: "~{memory_gb} GB"
    }
}
