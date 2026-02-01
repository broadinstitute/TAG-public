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
        # Main Outputs
        File viral_bam = FilterViralBam.out_bam
        File viral_bai = FilterViralBam.out_bai
        
        # Aligned Read Counts as Integers
        Int input_aligned_read_count = FilterViralBam.input_aligned_count
        Int viral_aligned_read_count = FilterViralBam.viral_aligned_count
        
        # Diagnostic Outputs
        File viral_names = FilterViralBam.viral_names
        File viral_fastq = FilterViralBam.viral_fastq
        File bbduk_stats = FilterViralBam.bbduk_stats
    }
}

task FilterViralBam {
    input {
        File bam_file
        File reference_fasta
        String basename
        
        # Resource Configuration
        Int threads = 8
        Int memory_gb = 16
        Int preemptible_attempts

        # Disk size calculation: 20x to handle large BAM->FASTA expansion
        Int disk_size_gb = ceil(20 * size(bam_file, "GB")) + 20

        String docker_image = "fleharty/viral-bam-filter:v2" 
    }

    command <<<
        set -e
        
        echo "--- DIAGNOSTIC START ---"
        echo "Processing Sample: ~{basename}"
        
        # -------------------------------------------------------
        # STEP 0: Count Aligned Reads in Input BAM
        # -------------------------------------------------------
        echo "Step 0: Counting aligned reads in input..."
        # -c counts, -F 4 excludes unmapped reads
        samtools view -c -F 4 "~{bam_file}" > input_aligned.txt
        echo "Input Aligned Count: $(cat input_aligned.txt)"

        # -------------------------------------------------------
        # STEP 1: Convert BAM to FASTA
        # -------------------------------------------------------
        echo "Step 1: Converting BAM to FASTA..."
        date
        samtools fastq -@ ~{threads} "~{bam_file}" > "~{basename}.fasta"
        
        echo "Step 1 Finished."
        date

        # -------------------------------------------------------
        # STEP 2: Run BBDuk (Viral Filtering)
        # -------------------------------------------------------
        echo "Step 2: Running BBDuk..."
        
        bbduk.sh -Xmx~{memory_gb - 4}g \
            in="~{basename}.fasta" \
            ref="~{reference_fasta}" \
            outm="~{basename}.viral.fastq" \
            k=19 hdist=0 rcomp=t \
            threads=~{threads} tbo=f tpe=f 2> bbduk_stats.txt

        echo "DIAGNOSTIC: BBDuk Finished. Stats:"
        cat bbduk_stats.txt
        
        # Cleanup intermediate FASTA to save disk space
        rm "~{basename}.fasta"
        
        # -------------------------------------------------------
        # STEP 3: Extract Read Names
        # -------------------------------------------------------
        echo "Step 3: Extracting read names..."
        
        # This strips /1 or /2 before saving the names to ensure match in original BAM
        grep '^@' "~{basename}.viral.fastq" | sed 's/^@//' | sed 's/\/[12]$//' | cut -d ' ' -f1 | sort -u > "~{basename}.viral.names"

        echo "DIAGNOSTIC: Count of Unique Viral Read Names Found: $(wc -l < "~{basename}.viral.names")"
        
        # -------------------------------------------------------
        # STEP 4: Subset Original BAM
        # -------------------------------------------------------
        echo "Step 4: Subsetting original BAM..."
        
        # -N filters the BAM based on the list of names
        samtools view -N "~{basename}.viral.names" -b "~{bam_file}" > "~{basename}.viral.bam"
        samtools index "~{basename}.viral.bam"
        
        # -------------------------------------------------------
        # STEP 5: Count Aligned Reads in Output Viral BAM
        # -------------------------------------------------------
        echo "Step 5: Counting aligned reads in output..."
        samtools view -c -F 4 "~{basename}.viral.bam" > viral_aligned.txt
        echo "Viral Aligned Count: $(cat viral_aligned.txt)"

        echo "Done."
    >>>

    output {
        File out_bam = "~{basename}.viral.bam"
        File out_bai = "~{basename}.viral.bam.bai"
        
        # Capturing counts as Integers for the workflow output
        Int input_aligned_count = read_int("input_aligned.txt")
        Int viral_aligned_count = read_int("viral_aligned.txt")

        # Debugging / Diagnostic Outputs
        File viral_names = "~{basename}.viral.names"
        File viral_fastq = "~{basename}.viral.fastq"
        File bbduk_stats = "bbduk_stats.txt"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: "~{memory_gb} GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: preemptible_attempts
    }
}
