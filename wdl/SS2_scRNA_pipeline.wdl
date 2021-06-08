import "SmartSeq_wdls/hisat2_descriptor.wdl" as run_hisat2
import "SmartSeq_wdls/hisat2_rsem_descriptor.wdl" as run_hisat2_rsem

## main pipeline:
## QC track: HISAT2+Picard
## this track will produce aligned bam file and a set of QC metrics
## rsem quantification track: HISAT2+RSEM
## this track involves hisat2 transcriptome alignment and RSEM gene expression estimation.

## Snapshot 23

workflow SmartSeq2SingleCell {

  # load annotation
  File gtf_file
  File genome_ref_fasta
  File rrna_intervals
  File gene_ref_flat
  #load index
  File hisat2_ref_index
  File hisat2_ref_trans_index
  File rsem_ref_index
  # ref index name
  String hisat2_ref_name
  String hisat2_ref_trans_name
  # samples
  String stranded
  String sample_name
  String output_name
  String smid

  String? ss2_docker
  String? ss2_adapter_qc_docker

  File fastq1
  File fastq2
  # adapter task
  Boolean check_adapter
  

  call run_hisat2.RunHisat2Pipeline as qc {
    input:
      fastq_read1 = fastq1,
      fastq_read2 = fastq2, 
      gtf = gtf_file,
      stranded = stranded,
      ref_fasta = genome_ref_fasta,
      rrna_interval = rrna_intervals,
      ref_flat = gene_ref_flat,
      hisat2_ref = hisat2_ref_index,
      hisat2_ref_name = hisat2_ref_name,
      sample_name = sample_name,
      output_prefix = output_name + "_qc"
      
  }

  call run_hisat2_rsem.RunHisat2RsemPipeline as data {
    input:
      fastq_read1 = fastq1, 
      fastq_read2 = fastq2, 
      hisat2_ref_trans = hisat2_ref_trans_index,
      hisat2_ref_trans_name = hisat2_ref_trans_name,
      rsem_genome = rsem_ref_index,
      output_prefix = output_name + "_rsem",
      sample_name = sample_name
  }
  if(check_adapter){
    call AdapterQC as adapter {
      input:
          ss2_adapter_qc_docker = select_first([ss2_adapter_qc_docker, "us.gcr.io/tag-team-160914/tag-tools:1.0.0"]),
          fastq1 = fastq1, 
          fastq2 = fastq2, 
          output_prefix = output_name
    }
  }

    call ExtractQC_metrics {
      input:
          ss2_docker = select_first([ss2_docker, "us.gcr.io/tag-team-160914/smartseq2_metrics:v1"]),
          rsem_gene_results = data.rsem_gene_results,
          dedup_metrics = qc.dedup_metrics,
          alignment_summary_metrics = qc.alignment_summary_metrics,
          rna_metrics = qc.rna_metrics,
          output_prefix = output_name

    }

  output {
    ## qc outputs
    File aligned_bam = qc.aligned_bam
    File? alignment_summary_metrics = qc.alignment_summary_metrics
    File? bait_bias_detail_metrics = qc.bait_bias_detail_metrics
    File? bait_bias_summary_metrics = qc.bait_bias_summary_metrics
    File? base_call_dist_metrics = qc.base_call_dist_metrics
    File? base_call_pdf = qc.base_call_pdf
    File? dedup_metrics = qc.dedup_metrics
    File? error_summary_metrics = qc.error_summary_metrics
    File? gc_bias_detail_metrics = qc.gc_bias_detail_metrics
    File? gc_bias_dist_pdf = qc.gc_bias_dist_pdf
    File? gc_bias_summary_metrics = qc.gc_bias_summary_metrics
    File? insert_size_hist = qc.insert_size_hist
    File? insert_size_metrics = qc.insert_size_metrics
    File hisat2_logfile = qc.logfile
    File hisat2_metfile = qc.metfile
    File? pre_adapter_details_metrics = qc.pre_adapter_details_metrics
    File? quality_by_cycle_metrics = qc.quality_by_cycle_metrics
    File? quality_by_cycle_pdf = qc.quality_by_cycle_pdf
    File? quality_distribution_dist_pdf = qc.quality_distribution_dist_pdf
    File? quality_distribution_metrics = qc.quality_distribution_metrics
    File? rna_coverage = qc.rna_coverage
    File? rna_metrics = qc.rna_metrics
    ## data outputs
    File aligned_trans_bam = data.aligned_trans_bam
    File hisat2tran_logfile = data.logfile
    File hisat2tran_metfile = data.metfile
    File? rsem_cnt_log = data.rsem_cnt_log
    File? rsem_gene_results = data.rsem_gene_results
    File? rsem_isoform_results = data.rsem_isoform_results
    File? rsem_model_log = data.rsem_model_log
    File? rsem_theta_log = data.rsem_theta_log
    File? rsem_time_log = data.rsem_time_log
    ## adapter outputs
    File? adapter_content_metrics = adapter.adapter_content_metrics
    File? read1_adapter_html = adapter.read1_fastqc_html
    File? read2_adapter_html = adapter.read2_fastqc_html

    #ExtractQC_metrics outputs
    Float pct_duplication = ExtractQC_metrics.merged_metrics["PERCENT_DUPLICATION"]
    Int estimated_library_size = ceil(ExtractQC_metrics.merged_metrics["ESTIMATED_LIBRARY_SIZE"])
    Int total_reads = ceil(ExtractQC_metrics.merged_metrics["TOTAL_READS"])
    Float pct_aligned = ExtractQC_metrics.merged_metrics["PCT_PF_READS_ALIGNED_PAIR"]
    Float pct_aligned_read1 = ExtractQC_metrics.merged_metrics["PCT_ALIGNED_READ1"]
    Float pct_aligned_read2 = ExtractQC_metrics.merged_metrics["PCT_ALIGNED_READ2"]
    Float pct_ribo = ExtractQC_metrics.merged_metrics["PCT_RIBOSOMAL_BASES"]
    Float pct_coding = ExtractQC_metrics.merged_metrics["PCT_CODING_BASES"]
    Float pct_utr = ExtractQC_metrics.merged_metrics["PCT_UTR_BASES"]
    Float pct_intronic = ExtractQC_metrics.merged_metrics["PCT_INTRONIC_BASES"]
    Float pct_intergenic = ExtractQC_metrics.merged_metrics["PCT_INTERGENIC_BASES"]
    Float pct_mRNA = ExtractQC_metrics.merged_metrics["PCT_MRNA_BASES"]
    Float pct_usable = ExtractQC_metrics.merged_metrics["PCT_USABLE_BASES"]
    Float median_cv_coverage = ExtractQC_metrics.merged_metrics["MEDIAN_CV_COVERAGE"]
    Float median_5prime_bias = ExtractQC_metrics.merged_metrics["MEDIAN_5PRIME_BIAS"]
    Float median_3prime_bias = ExtractQC_metrics.merged_metrics["MEDIAN_3PRIME_BIAS"]
    Int genes_detected = ceil(ExtractQC_metrics.merged_metrics["GENES_DETECTED"])
    Float pct_mito = ExtractQC_metrics.merged_metrics["PCT_MITO"]

  }
}

task AdapterQC {
  File fastq1
    File fastq2
    File adapter_script
    String output_prefix
    
    String ss2_adapter_qc_docker

    Float? mem = 4
    
    Float fastq_size = size(fastq1, "GB") + size(fastq2, "GB")
    Int? increase_disk_size = 50
    Int disk_size = ceil(fastq_size * 10)  + increase_disk_size

    String basename1 = basename(fastq1,".fastq.gz")
    String basename2 = basename(fastq2,".fastq.gz")

    command <<<
      perl /usr/tag/scripts/FastQC/fastqc --extract ${fastq1} ${fastq2} -o ./

        python ${adapter_script} ${basename1} ${basename2} ${output_prefix}
    >>>
    runtime {
      docker: ss2_adapter_qc_docker
      disks: "local-disk " + disk_size + " HDD"
      memory: mem + "GB"
      cpu: "1"
    }
    output {
      File adapter_content_metrics = "${output_prefix}.adapter_content_metrics.txt"
        File read1_fastqc_html = "${basename1}_fastqc.html"
        File read2_fastqc_html = "${basename2}_fastqc.html"
    }
}

task ExtractQC_metrics {
    File rsem_gene_results
    File dedup_metrics
    File alignment_summary_metrics
    File rna_metrics

    String ss2_docker

    String output_prefix

    Float? mem = 4
    

    command <<<

    paste <(cat ${dedup_metrics} | sed '/^#/d' | cut -f 9,10 | grep -v ^$) \
    <(cat ${alignment_summary_metrics} | sed '/^#/d' | cut -f 2,7 | sed -n -e 2p -e 5p | sed 's/PCT_PF_READS_ALIGNED/PCT_PF_READS_ALIGNED_PAIR/g') \
    <(cat ${alignment_summary_metrics} | sed '/^#/d' | cut -f 2,7 | sed -n -e 2p -e 3p | cut -f 2 | sed 's/PCT_PF_READS_ALIGNED/PCT_ALIGNED_READ1/g') \
    <(cat ${alignment_summary_metrics} | sed '/^#/d' | cut -f 2,7 | cut -f 2 | sed -n -e 2p -e 4p | sed 's/PCT_PF_READS_ALIGNED/PCT_ALIGNED_READ2/g') \
    <(cat ${rna_metrics} | sed '/^#/d' | cut -f 16-22,24-26 | grep -v ^$) \
    <(python3 /scripts/rsem_gene_results_metrics.py --rsem_genes_results ${rsem_gene_results}) > metrics_in_columns.txt
    grep ^[A-Z] ${output_prefix}.metric_in_rows.txt | awk 'BEGIN{FS="\t"}{for(i=1; i<=NF; i++){a[NR,i] = $i}} NF>p { p = NF }END{for(j=1; j<=p; j++){str=a[1,j]; for(i=2; i<=NR; i++){if(a[i,j] == ""){a[i,j]=-1} str=str"\t"a[i,j]};print str}}' metrics_in_columns.txt > ${output_prefix}.metric_in_rows.txt

    >>>
    runtime {
      docker: ss2_docker
      disks: "local-disk 25 HDD"
      memory: mem + "GB"
      cpu: "1"
    }

    output {

        File merged_metrics_file = "${output_prefix}.metric_in_rows.txt"
        Map[String, Float] merged_metrics = read_map("${output_prefix}.metric_in_rows.txt")

    }


}
