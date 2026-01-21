version 1.0

workflow SplitVcfWorkflow {
   call SplitVCFs
}

task SplitVCFs {

   input {
      String splitvcf_docker = "us.gcr.io/broad-dsde-methods/liquidbiopsy:0.0.4.3"
      String basename
      File? gatk_override
      File vcf_file
      File vcf_idx
      File reference_fasta
      File reference_fasta_idx
      File reference_dict
      Int? preemptible
      Int? memory
      Int? disk_pad
   }

   Int ref_size = ceil(size(reference_fasta, "GB") + size(reference_fasta_idx, "GB") + size(reference_dict, "GB"))
   Int disk_size = ceil(size(vcf_file, "GB") * 2) + ceil(size(vcf_idx, "GB")) + ref_size + select_first([disk_pad, 0])
   Int mem = select_first([memory, 16])

   command <<<

      set -e
      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      # Include all snp and mnp calls in the snp vcf
      gatk --java-options "-Xmx15G" SelectVariants \
         -V ~{vcf_file} \
         -O "~{basename}.snp.vcf.gz" \
         -R ~{reference_fasta} \
         -select-type SNP \
         -select-type MNP

      #passing only SNPs
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "~{basename}.snp.vcf.gz" \
         -O "~{basename}.snp.passing.vcf.gz" \
         -R ~{reference_fasta} \
         --exclude-filtered \
         -select-type SNP

      passing_SNP="$(gatk CountVariants -V ~{basename}.snp.passing.vcf.gz | tail -1)"
      echo "$passing_SNP" > passing_snp.txt

      #passing only MNPs
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "~{basename}.snp.vcf.gz" \
         -O "~{basename}.mnp.passing.vcf.gz" \
         -R ~{reference_fasta} \
         --exclude-filtered \
         -select-type MNP

      passing_MNP="$(gatk CountVariants -V ~{basename}.mnp.passing.vcf.gz | tail -1)"
      echo "$passing_MNP" > passing_mnp.txt

      #filtered only SNPs
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "~{basename}.snp.vcf.gz" \
         -XL "~{basename}.snp.passing.vcf.gz" \
         -O "~{basename}.snp.filtered.vcf.gz" \
         -select-type SNP \
         -R ~{reference_fasta}

      filtered_SNP="$(gatk CountVariants -V ~{basename}.snp.filtered.vcf.gz | tail -1)"
      echo "$filtered_SNP" > filtered_snp.txt

      #filtered only MNPs
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "~{basename}.snp.vcf.gz" \
         -XL "~{basename}.mnp.passing.vcf.gz" \
         -O "~{basename}.mnp.filtered.vcf.gz" \
         -select-type MNP \
         -R ~{reference_fasta}

      filtered_MNP="$(gatk CountVariants -V ~{basename}.mnp.filtered.vcf.gz | tail -1)"
      echo "$filtered_MNP" > filtered_mnp.txt

      gatk --java-options "-Xmx15G" SelectVariants \
         -V ~{vcf_file} \
         -O "~{basename}.indel.vcf.gz" \
         -R ~{reference_fasta} \
         -select-type INDEL \
         -select-type MIXED

      #passing only
      gatk --java-options "-Xmx15G" SelectVariants \
         -V "~{basename}.indel.vcf.gz" \
         -O "~{basename}.indel.passing.vcf.gz" \
         -R ~{reference_fasta} \
         --exclude-filtered

      passing_INDEL="$(gatk CountVariants -V ~{basename}.indel.passing.vcf.gz | tail -1)"
      echo "$passing_INDEL" > passing_indel.txt

      #filtered only INDELs
      bcftools isec -C -O v -o ~{basename}.indel.filtered.vcf -w1 ~{basename}.indel.vcf.gz ~{basename}.indel.passing.vcf.gz

      filtered_INDEL="$(gatk CountVariants -V ~{basename}.indel.filtered.vcf | tail -1)"
      echo "$filtered_INDEL" > filtered_indel.txt

   >>>

   runtime {
      docker: splitvcf_docker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + " GB"
      preemptible: select_first([preemptible, 3])
   }
   output {

      File output_snp_vcf = "${basename}.snp.vcf.gz"
      File output_snp_vcf_idx = "${basename}.snp.vcf.gz.tbi"

      Int passing_SNP = read_int("passing_snp.txt")
      Int filtered_SNP = read_int("filtered_snp.txt")

      Int passing_MNP = read_int("passing_mnp.txt")
      Int filtered_MNP = read_int("filtered_mnp.txt")

      File output_indel_vcf = "${basename}.indel.vcf.gz"
      File output_indel_vcf_idx = "${basename}.indel.vcf.gz.tbi"

      Int passing_INDEL = read_int("passing_indel.txt")
      Int filtered_INDEL = read_int("filtered_indel.txt")

   }
}