version 1.0

workflow FilterViralReads {
    input {
        File input_bam
        File viral_reference
        String sample_name
        Int preemptible_attempts = 1
    }

    call FilterViralBam {
        input:
            bam_file = input_bam,
            reference_fasta = viral_reference,
            basename = sample_name,
            preemptible_attempts = preemptible_attempts
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
        Int preemptible_attempts

        # Calculate disk: (20 * size of BAM in GB) + 10GB buffer for tools/logs
        # The reason for 20x size is that it converts the bam to a fasta file.
        Int disk_size_gb = ceil(10 * size(bam_file, "GB")) + 10

        # Replace this string with the tag of the image you built in Part 1
        String docker_image = "fleharty/viral-bam-filter:v2" 
    }

    command <<<
        set -e
        echo "df -h"
        df -h 
        
        echo "df -h ."
        df -h .

        echo "pwd"
        pwd

        echo "du -h"
        du -h 

        echo "Processing Sample: ~{basename}"
        
        # 1. Convert BAM to FASTA
        # We use ~{threads} to utilize the requested CPU cores
        echo "Step 1: Convert Bam to FASTA"
        samtools fastq -@ ~{threads} "~{bam_file}" | pv -i 10 2> samtools_progress.log > "~{basename}.fasta" &

        # 2. Run BBDuk
        echo "Step 2: Run BBDuk"
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
        echo "Step 3: Extract read names"
        grep '^@' "~{basename}.viral.fastq" | sed 's/^@//' | cut -d ' ' -f1 | sort -u > "~{basename}.viral.names"

        # 4. Subset original BAM using the read names
        echo "Step 4: Subset original BAM using the read names"
        samtools view -N "~{basename}.viral.names" -b "~{bam_file}" > "~{basename}.viral.bam"
        samtools index "~{basename}.viral.bam"

        # Final Cleanup
        echo "Step 5: Final Cleanup"
        rm "~{basename}.viral.names"
        rm "~{basename}.viral.fastq"

        echo "Done"
    >>>

    output {
        File out_bam = "~{basename}.viral.bam"
        File out_bai = "~{basename}.viral.bam.bai"
        File progress_log = "samtools_progress.log"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: "~{memory_gb} GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: preemptible_attempts
    }
}
