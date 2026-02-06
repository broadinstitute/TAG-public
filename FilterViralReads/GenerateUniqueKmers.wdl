version 1.0

workflow GenerateUniqueKmers {
    input {
        File viral_vector_fasta
        File host_reference_fasta
        String output_basename
        Int preemptible_attempts = 1
    }

    meta {
        description: "Generates a reference set of viral k-mers that do not exist in the host genome using BBDuk subtraction."
        author: "fleharty"
        email: "fleharty@broadinstitute.org"
    }

    call SubtractHostKmers {
        input:
            viral_vector = viral_vector_fasta,
            host_reference = host_reference_fasta,
            basename = output_basename,
            preemptible_attempts = preemptible_attempts
    }

    output {
        File unique_viral_kmers = SubtractHostKmers.out_unique_fasta
        File matched_host_kmers = SubtractHostKmers.out_matched_fasta
    }
}

task SubtractHostKmers {
    input {
        File viral_vector
        File host_reference
        String basename
        
        # Resource Configuration
        Int threads = 16
        Int memory_gb = 40  # Higher memory for hosting human genome index
        Int preemptible_attempts
        
        # BBDuk Parameters
        Int kmer_size = 19
        Int hdist = 0

        # Disk size: Host reference can be large (~3GB) + viral files
        Int disk_size_gb = ceil(size(host_reference, "GB") + size(viral_vector, "GB")) + 50

        String docker_image = "fleharty/viral-bam-filter:v2" 
    }

    command <<<
        set -e
        
        echo "--- STEP 1: K-merizing the Viral Vector ---"
        bbduk.sh -Xmx8g \
            in="~{viral_vector}" \
            out="viral_k~{kmer_size}.fasta" \
            k=~{kmer_size} \
            mink=~{kmer_size}

        echo "--- STEP 2: Subtracting Host K-mers ---"
        # Using memory_gb - 4 for the Java heap to leave overhead for the OS
        bbduk.sh -Xmx~{memory_gb - 4}g \
            in="viral_k~{kmer_size}.fasta" \
            ref="~{host_reference}" \
            outu="~{basename}_unique.fasta" \
            outm="~{basename}_matched_host.fasta" \
            k=~{kmer_size} \
            hdist=~{hdist} \
            rcomp=t \
            trd=t \
            threads=~{threads}
            
        echo "--- Generation Complete ---"
    >>>

    output {
        File out_unique_fasta = "~{basename}_unique.fasta"
        File out_matched_fasta = "~{basename}_matched_host.fasta"
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: "~{memory_gb} GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: preemptible_attempts
    }
}
