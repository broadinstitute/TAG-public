version 1.0

workflow FilterViralReads {
    input {
        File input_bam
        File viral_reference
        String sample_name
        Int preemptible_attempts = 1
    }

    meta {
        description: "Extracts viral-specific reads from a BAM file using BBDuk k-mer matching. Identifies reads (and mates) matching the provided viral reference and outputs a subsetted BAM." [cite: 12, 14]
        author: "fleharty"
        email: "fleharty@broadinstitute.org"
    }

    parameter_meta {
        input_bam: "The original BAM file to be filtered; can be aligned to any host reference." [cite: 7]
        viral_reference: "FASTA file containing viral sequences. K-mers should be unique to the virus to avoid host contamination." 
        sample_name: "Prefix used for all output files." [cite: 3]
    }

    call FilterViralBam {
        input:
            bam_file = input_bam,
            reference_fasta = viral_reference,
            basename = sample_name, [cite: 3]
            preemptible_attempts = preemptible_attempts
    }
    
    output {
        File viral_bam = FilterViralBam.out_bam
        File viral_bai = FilterViralBam.out_bai
        Int input_aligned_read_count = FilterViralBam.input_aligned_count [cite: 4]
        Int viral_aligned_read_count = FilterViralBam.viral_aligned_count [cite: 4]
        File viral_names = FilterViralBam.viral_names
        File viral_fastq = FilterViralBam.viral_fastq
        File bbduk_stats = FilterViralBam.bbduk_stats
    }
}

task FilterViralBam {
    input {
        File bam_file
        File reference_fasta
        String basename [cite: 5]
        
        # Resource Configuration
        Int threads = 8
        Int memory_gb = 16
        Int preemptible_attempts

        # Configurable BBDuk parameters
        Int kmer_size = 19 
        Int hdist = 0 

        # Disk size calculation: Adjusted to 10x for safety without being excessive
        Int disk_size_gb = ceil(10 * size(bam_file, "GB")) + 20 [cite: 5]

        String docker_image = "fleharty/viral-bam-filter:v2" 
    }

    parameter_meta {
        threads: "Number of CPUs to allocate for samtools and BBDuk." [cite: 8, 10]
        memory_gb: "Total memory; note that BBDuk is allocated memory_gb - 4GB." 
        kmer_size: "K-mer length for matching (default 19)." 
        hdist: "Hamming distance for k-mer matching; 0 requires exact matches." 
    }

    command <<<
        set -e
        
        echo "Step 0: Counting aligned reads in input..."
        samtools view -c -F 4 "~{bam_file}" > input_aligned.txt [cite: 7]

        echo "Step 1: Converting BAM to FASTQ..."
        samtools fastq -@ ~{threads} "~{bam_file}" > "~{basename}.fastq" [cite: 8]
        
        echo "Step 2: Running BBDuk..."
        # BBDuk uses k-mer matching to identify viral reads 
        bbduk.sh -Xmx~{memory_gb - 4}g \
            in="~{basename}.fastq" \
            ref="~{reference_fasta}" \
            outm="~{basename}.viral.fastq" \
            k=~{kmer_size} hdist=~{hdist} rcomp=t \
            threads=~{threads} tbo=f tpe=f 2> bbduk_stats.txt 
        
        rm "~{basename}.fastq"
        
        echo "Step 3: Extracting unique read names..."
        # Strips /1 and /2 to ensure both mates are captured from the original BAM [cite: 12]
        grep '^@' "~{basename}.viral.fastq" | sed 's/^@//' | sed 's/\/[12]$//' | cut -d ' ' -f1 | sort -u > "~{basename}.viral.names" [cite: 12, 13]

        echo "Step 4: Subsetting original BAM..."
        # -N uses the name list to pull full read records [cite: 14]
        samtools view -N "~{basename}.viral.names" -b "~{bam_file}" > "~{basename}.viral.bam" [cite: 14]
        samtools index "~{basename}.viral.bam"
        
        echo "Step 5: Counting aligned reads in output..."
        samtools view -c -F 4 "~{basename}.viral.bam" > viral_aligned.txt [cite: 15]
    >>>

    output {
        File out_bam = "~{basename}.viral.bam"
        File out_bai = "~{basename}.viral.bam.bai"
        Int input_aligned_count = read_int("input_aligned.txt") [cite: 16]
        Int viral_aligned_count = read_int("viral_aligned.txt") [cite: 16]
        File viral_names = "~{basename}.viral.names"
        File viral_fastq = "~{basename}.viral.fastq" [cite: 17]
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
