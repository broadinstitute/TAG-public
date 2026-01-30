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
        
        Int threads = 8
        Int memory_gb = 16
        Int preemptible_attempts

        # Disk size calculation
        Int disk_size_gb = ceil(10 * size(bam_file, "GB")) + 10

        String docker_image = "fleharty/viral-bam-filter:v2" 
    }

    command <<<
        set -e
        
        echo "--- DIAGNOSTIC: STARTING ---"
        echo "Input BAM Size:"
        ls -lh "~{bam_file}"

        # 1. Convert BAM to FASTA
        # REMOVED THE '&' to fix race condition
        echo "Step 1: Converting BAM to FASTA..."
        samtools fastq -@ ~{threads} "~{bam_file}" > "~{basename}.fasta"
        
        echo "DIAGNOSTIC: FASTA generated. Size:"
        ls -lh "~{basename}.fasta"
        echo "DIAGNOSTIC: First 5 lines of FASTA:"
        head -n 5 "~{basename}.fasta"

        # 2. Run BBDuk
        echo "Step 2: Running BBDuk..."
        # Capturing BBDuk stats to a file so we can see them in outputs
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
        
        # Check if BBDuk found 0 reads
        if [ ! -s "~{basename}.viral.fastq" ]; then
            echo "WARNING: BBDuk output is empty! No viral reads found."
        fi

        # 3. Extract read names
        echo "Step 3: Extracting read names..."
        # We save this output to verify if suffixes like '/1' are contaminating the names
        grep '^@' "~{basename}.viral.fastq" | sed 's/^@//' | cut -d ' ' -f1 | sort -u > "~{basename}.viral.names"

        echo "DIAGNOSTIC: Read Names File Size:"
        ls -lh "~{basename}.viral.names"
        echo "DIAGNOSTIC: Count of Viral Reads Found:"
        wc -l "~{basename}.viral.names"
        echo "DIAGNOSTIC: First 5 Read Names (Check for /1 or /2 suffixes!):"
        head -n 5 "~{basename}.viral.names"

        # 4. Subset original BAM
        echo "Step 4: Subsetting BAM..."
        samtools view -N "~{basename}.viral.names" -b "~{bam_file}" > "~{basename}.viral.bam"
        samtools index "~{basename}.viral.bam"
        
        echo "DIAGNOSTIC: Final BAM Size:"
        ls -lh "~{basename}.viral.bam"
        
        # Check if final BAM has reads
        echo "DIAGNOSTIC: Checking final BAM read count:"
        samtools view -c "~{basename}.viral.bam"

        # Cleanup (Only delete the massive FASTA, keep others for debugging)
        rm "~{basename}.fasta"
    >>>

    output {
        File out_bam = "~{basename}.viral.bam"
        File out_bai = "~{basename}.viral.bam.bai"
        
        # NEW DEBUGGING OUTPUTS
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
