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
        echo "Input BAM Size:"
        ls -lh "~{bam_file}"

        # -------------------------------------------------------
        # STEP 1: Convert BAM to FASTA
        # -------------------------------------------------------
        echo "Step 1: Converting BAM to FASTA..."
        date # Print timestamp so you can calculate duration later

        # Standard conversion. Silent until finished.
        samtools fastq -@ ~{threads} "~{bam_file}" > "~{basename}.fasta"
        
        echo "Step 1 Finished."
        date
        echo "DIAGNOSTIC: FASTA generated. Size:"
        ls -lh "~{basename}.fasta"
        echo "DIAGNOSTIC: First 4 lines of FASTA:"
        head -n 4 "~{basename}.fasta"

        # -------------------------------------------------------
        # STEP 2: Run BBDuk (Viral Filtering)
        # -------------------------------------------------------
        echo "Step 2: Running BBDuk..."
        
        # Saving BBDuk statistics to a text file for output
        bbduk.sh -Xmx~{memory_gb - 4}g \
            in="~{basename}.fasta" \
            ref="~{reference_fasta}" \
            outm="~{basename}.viral.fastq" \
            k=19 hdist=0 rcomp=t \
            threads=~{threads} tbo=f tpe=f 2> bbduk_stats.txt

        echo "DIAGNOSTIC: BBDuk Finished. Stats:"
        cat bbduk_stats.txt
        
        echo "DIAGNOSTIC: Viral FASTQ Size:"
        ls -lh "~{basename}.viral.fastq"
        
        # -------------------------------------------------------
        # STEP 3: Extract Read Names
        # -------------------------------------------------------
        echo "Step 3: Extracting read names..."
        
        # This strips /1 or /2 before saving the names
        grep '^@' "~{basename}.viral.fastq" | sed 's/^@//' | sed 's/\/[12]$//' | cut -d ' ' -f1 | sort -u > "~{basename}.viral.names"

        echo "DIAGNOSTIC: Read Names File Size:"
        ls -lh "~{basename}.viral.names"
        echo "DIAGNOSTIC: Count of Viral Reads Found:"
        wc -l "~{basename}.viral.names"
        
        echo "DIAGNOSTIC: First 10 Read Names (Check for /1 or /2 suffixes!):"
        head -n 10 "~{basename}.viral.names"

        # -------------------------------------------------------
        # STEP 4: Subset Original BAM
        # -------------------------------------------------------
        echo "Step 4: Subsetting original BAM..."
        
        # -N filters the BAM based on the list of names we just generated
        samtools view -N "~{basename}.viral.names" -b "~{bam_file}" > "~{basename}.viral.bam"
        samtools index "~{basename}.viral.bam"
        
        echo "DIAGNOSTIC: Final BAM Size:"
        ls -lh "~{basename}.viral.bam"
        
        echo "DIAGNOSTIC: Checking final BAM read count:"
        samtools view -c "~{basename}.viral.bam"

        # -------------------------------------------------------
        # Cleanup
        # -------------------------------------------------------
        echo "Cleaning up large FASTA file..."
        rm "~{basename}.fasta"
        echo "Done."
    >>>

    output {
        File out_bam = "~{basename}.viral.bam"
        File out_bai = "~{basename}.viral.bam.bai"
        
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
