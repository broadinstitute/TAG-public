## quantification pipeline
## hisat2 as aligner to align reads to transcriptome
## output: transcriptome aligned bam files
## rsem will estimate the expression levels for both gene/isoform
## rsem will estimate abundance based on maximum likelihood (MLE) and posterior mean estimate (PME)
## output: genes.results and isoform.results 

workflow RunHisat2RsemPipeline {
  File fastq_read1
  File fastq_read2
  File hisat2_ref_trans
  File rsem_genome
  String output_prefix
  String hisat2_ref_trans_name
  String sample_name
  ##variable to estimate disk size
  ## variables to estimate disk size
  Float hisat2_ref_size = size(hisat2_ref_trans, "GB")
  Float fastq_size = size(fastq_read1, "GB") +size(fastq_read2, "GB")
  Float rsem_ref_size = size(rsem_genome, "GB")
  Float bam_disk_multiplier = 10.0
  Int? increase_disk_size
  Int? Rsem_mem 
  Int? Rsem_exp_mem
  Int Rsem_memory = select_first([Rsem_mem, 5])
  Int Rsem_exp_memory = select_first([Rsem_exp_mem, 4])
  Int additional_disk = select_first([increase_disk_size, 10])
  Boolean zip_fq1 = if basename(fastq_read1,".gz") == basename(fastq_read1) then true else false
  Boolean zip_fq2 = if basename(fastq_read2,".gz") == basename(fastq_read2) then true else false
  
  call HISAT2rsem as Hisat2Trans {
    input:
      hisat2_ref = hisat2_ref_trans,
      fq1 = fastq_read1,
      fq2 = fastq_read2,
      zip_fq1 = zip_fq1,
      zip_fq2 = zip_fq2,
      ref_name = hisat2_ref_trans_name,
      sample_name = sample_name,
      output_name = output_prefix,
      disk_size = ceil(fastq_size * bam_disk_multiplier + hisat2_ref_size + additional_disk * 5.0),
      Rsem_memory = Rsem_memory

    }
  #Boolean reads_aligned = Hisat2Trans.reads_aligned  
  Float bam_size = size(Hisat2Trans.output_bam, "GB")
  #if (reads_aligned) {
  call RsemExpression as Rsem {
    input:
      trans_aligned_bam = Hisat2Trans.output_bam,
      rsem_genome = rsem_genome,
      rsem_out = output_prefix,
      disk_size = ceil(fastq_size * bam_disk_multiplier+rsem_ref_size+additional_disk * 2.0),
      Rsem_exp_memory = Rsem_exp_memory
    }
  #}
  output {
    File aligned_trans_bam = Hisat2Trans.output_bam
    File metfile = Hisat2Trans.metfile
    File logfile = Hisat2Trans.logfile
    File? rsem_gene_results = Rsem.rsem_gene
    File? rsem_isoform_results = Rsem.rsem_isoform
    File? rsem_time_log = Rsem.rsem_time
    File? rsem_cnt_log = Rsem.rsem_cnt
    File? rsem_model_log = Rsem.rsem_model
    File? rsem_theta_log = Rsem.rsem_theta
  }
}


## run hisat2 for rsem
## increase gap alignment penalty to avoid gap alignment
## --mp 1,1 --np 1 --score-min L,0,-0.1 is default paramesters when rsem runs alignment by using bowtie2/Bowtie
## --mp 1,1 and --np 1 will reduce mismatching penalty to 1 for all.
## with no-splice-alignment no-softclip no-mixed options on, HISAT2 will only output concordant alignment without soft-cliping
## --rdg 99999999,99999999 and --rfg 99999999,99999999 will give an infinity penalty to alignment with indel.As results
## no indel/gaps in alignments 
task HISAT2rsem {
  File hisat2_ref
  File fq1
  File fq2
  Boolean zip_fq1
  Boolean zip_fq2
  String ref_name
  String output_name
  String sample_name
  Int disk_size
  Int Rsem_memory

  command {
    set -e

    # zip fastq files if necessary.
    if [ "${zip_fq1}" == "true" ]; then
        FQ1=${fq1}.fastq.gz
        gzip ${fq1} -c > ${fq1}.fastq.gz
    else
        FQ1=${fq1}
    fi
    echo $FQ1

    if [ "${zip_fq2}" == "true" ]; then
        FQ2=${fq2}.fastq.gz
        gzip ${fq2} -c > ${fq2}.fastq.gz
    else
        FQ2=${fq2}
    fi
    echo $FQ2

    tar -zxf "${hisat2_ref}"
    hisat2 -t \
      -x ${ref_name}/${ref_name} \
      -1 $FQ1 \
      -2 $FQ2 \
      --rg-id=${sample_name} --rg SM:${sample_name} --rg LB:${sample_name} \
      --rg PL:ILLUMINA --rg PU:${sample_name} \
      --new-summary --summary-file ${output_name}.log \
      --met-file ${output_name}.hisat2.met.txt --met 5 \
      -k 10 \
      --mp 1,1 \
      --np 1 \
      --score-min L,0,-0.1 \
      --secondary \
      --no-mixed \
      --no-softclip \
      --no-discordant \
      --rdg 99999999,99999999 \
      --rfg 99999999,99999999 \
      --no-spliced-alignment \
      --seed 12345 \
      -p 4 -S ${output_name}.sam 
    samtools view -bS  "${output_name}.sam" > "${output_name}.bam"
	
    # This section is to check for whether we're actually getting alignments.
    # We grep the log file, and if it has lines about the alignments, we can check that we got more than 0
    # We then write a boolean to stdout, which gets read into a Boolean variable in the outputs section. This boolean tells the wdl whether to run the next steps or not.
	#if [ `grep ' 1 time: ' "${output_name}.log" -c` -ne 0 ]; then
    #	    if [ `grep ' 1 time: ' "${output_name}.log" | head -n 1| cut -f 2 -d ':' | cut -f 1 -d '(' | tr -d [:space:]` -eq 0 ]; then
    #    	        echo "false"
    #        	    else echo "true"
    #    	fi
	#	else echo "false"

	#fi
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    memory: Rsem_memory + "GB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "4"
    preemptible: 5
  }
  output {
	#Boolean reads_aligned = read_boolean(stdout())
    File logfile = "${output_name}.log"
    File metfile = "${output_name}.hisat2.met.txt"
    File output_bam = "${output_name}.bam"
  }
}


## RSEM will estimate gene/isoform quantification
## --bam: input is bam file
## -p 4: run on multiple threads
## --time: report running time
## --seed: report deterministic results
## --calc-pme will do posterior mean estimation
## --single-cell-prior in pme estimation, use Dirichlet(0.1) as the prior 
## which encourage the sparsity of the expression levels
## of note, the --single-cell-prior has not been tested yet and there are
## more investigation or improvement requested.  
task RsemExpression {
  File trans_aligned_bam
  File rsem_genome
  String rsem_out
  Int disk_size
  Int Rsem_exp_memory
  command {
    set -e
  
    tar -xvf ${rsem_genome}
    rsem-calculate-expression \
      --bam \
      --paired-end \
       -p 4 \
      --time --seed 555 \
      --calc-pme \
      --single-cell-prior \
      ${trans_aligned_bam} \
      rsem/rsem_trans_index  \
      "${rsem_out}" 
  }
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-rsem:v0.2.2-1.3.0"
    memory: Rsem_exp_memory + "GB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "4"
    preemptible: 5
  }
  output {
    File rsem_gene = "${rsem_out}.genes.results"
    File rsem_isoform = "${rsem_out}.isoforms.results"
    File rsem_time = "${rsem_out}.time"
    File rsem_cnt = "${rsem_out}.stat/${rsem_out}.cnt"
    File rsem_model = "${rsem_out}.stat/${rsem_out}.model"
    File rsem_theta = "${rsem_out}.stat/${rsem_out}.theta"
  }
}
