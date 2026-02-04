version 1.0

workflow FilterViralReads {
    input {
        File input_bam
        File viral_reference
        String sample_name
        Int preemptible_attempts = 1
    }

    meta {
        description: "Extracts viral-specific reads from a BAM file using BBDuk k-mer matching. Identifies reads (and mates) matching the provided viral reference and outputs a subsetted BAM."
        author: "fleharty"
        email: "fleharty@broadinstitute.org"
    }

    parameter_meta {
        input_bam: "The original BAM file to be filtered; can be aligned to any host reference."
        viral_reference: "FASTA file containing viral sequences. K-mers should be unique to the virus to avoid host contamination."
        sample_name: "Prefix used for all output files."
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

        # Configurable BBDuk parameters
        Int kmer_size = 19
        Int hdist = 0

        # Disk size calculation: 10x to handle BAM->FASTQ expansion safely
        Int disk_size_gb = ceil(10 * size(bam_file, "GB")) + 20

        String docker_image = "fleharty/viral-bam-filter:v2" 
    }

    parameter_meta {
        threads: "Number of CPUs to allocate for samtools and BBDuk."
        memory_gb: "Total memory; note that BBDuk is allocated memory_gb - 4GB."
        kmer_size: "K-mer length for matching (default 19)."
        hdist: "Hamming distance for k-mer matching; 0 requires exact matches."
    }

    command <<<
        set -e
        
        echo "--- DIAGNOSTIC START ---"
        echo "Processing Sample: ~{basename}"
        
        # STEP 0: Count Aligned Reads in Input BAM
        echo "Step 0: Counting aligned reads in input..."
        samtools view -c -F 4 "~{bam_file}" > input_aligned.txt
        echo "Input Aligned Count: $(cat input_aligned.txt)"

        # STEP 1: Convert BAM to FASTQ
        echo "Step 1: Converting BAM to FASTQ..."
        samtools fastq -@ ~{threads} "~{bam_file}" > "~{basename}.fastq"
        
        # STEP 2: Run BBDuk (Viral Filtering)
        echo "Step 2: Running BBDuk..."
        bbduk.sh -Xmx~{memory_gb - 4}g \
            in="~{basename}.fastq" \
            ref="~{reference_fasta}" \
            outm="~{basename}.viral.fastq" \
            k=~{kmer_size} hdist=~{hdist} rcomp=t \
            threads=~{threads} tbo=f tpe=f 2> bbduk_stats.txt

        echo "DIAGNOSTIC: BBDuk Finished. Stats:"
        cat bbduk_stats.txt
        
        # Cleanup intermediate FASTQ to save disk space
        rm "~{basename}.fastq"
        
        # STEP 3: Extract Read Names
        echo "Step 3: Extracting unique read names..."
        # Strips /1 or /2 suffixes to ensure both mates are captured in the original BAM subset
        grep '^@' "~{basename}.viral.fastq" | sed 's/^@//' | sed 's/\/[12]$//' | cut -d ' ' -f1 | sort -u > "~{basename}.viral.names"

        echo "DIAGNOSTIC: Unique Viral Read Names Found: $(wc -l < "~{basename}.viral.names")"
        
        # STEP 4: Subset Original BAM
        echo "Step 4: Subsetting original BAM..."
        samtools view -N "~{basename}.viral.names" -b "~{bam_file}" > "~{basename}.viral.bam"
        samtools index "~{basename}.viral.bam"
        
        # STEP 5: Count Aligned Reads in Output Viral BAM
        echo "Step 5: Counting aligned reads in output..."
        samtools view -c -F 4 "~{basename}.viral.bam" > viral_aligned.txt
        echo "Viral Aligned Count: $(cat viral_aligned.txt)"

        echo "Done."
    >>>

    output {
        File out_bam = "~{basename}.viral.bam"
        File out_bai = "~{basename}.viral.bam.bai"
        Int input_aligned_count = read_int("input_aligned.txt")
        Int viral_aligned_count = read_int("viral_aligned.txt")
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
