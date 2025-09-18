version 1.0

workflow SingleSampleCODEC_targeted {
    input {
        String sample_id
        File fastq1
        File fastq2
        File reference_fasta
        File reference_fasta_index
        File reference_dict
        File reference_pac
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_sa
        File? germline_bam
        File? germline_bam_index
        Int num_parallel
        String sort_memory       
        File eval_genome_interval = "gs://gptag/CODEC/GRCh38_notinalldifficultregions.interval_list"
        File eval_genome_bed = "gs://gptag/CODEC/GRCh38_notinalldifficultregions.bed"
        File? bait_intervals
        File? target_intervals
        Boolean is_captured_data
        Boolean create_vcf
    }
        call SplitFastq1 {
            input: 
                fastq_read1 = fastq1,
                nsplit = num_parallel,
                sample_id = sample_id
        }
        call SplitFastq2 {
            input: 
                fastq_read2 = fastq2,
                nsplit = num_parallel,
                sample_id = sample_id
        }
    scatter (split_index in range(num_parallel)) {
        String output_prefix = "${sample_id}_split.${split_index}"
        call Trim {
            input:
                read1 = SplitFastq1.split_read1[split_index],
                read2 = SplitFastq2.split_read2[split_index],
                output_prefix = output_prefix,
                split = split_index,
                sample_id = sample_id
        }
        call AlignRawTrimmed {
            input:
                bam_input = Trim.trimmed_bam,
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                reference_dict = reference_dict,
                reference_pac = reference_pac,
                reference_amb = reference_amb,
                reference_ann = reference_ann,
                reference_bwt = reference_bwt,
                reference_sa = reference_sa,
                sample_id = sample_id,
                split = split_index              
        }
        call ZipperBamAlignment {            
            input:
                mapped_bam = AlignRawTrimmed.aligned_bam,
                unmapped_bam = Trim.trimmed_bam,
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                reference_dict = reference_dict,
                sample_id = sample_id,
                split = split_index,
                sort_memory = sort_memory
            }
        }
        call MergeSplit {
            input:
                bam_files = ZipperBamAlignment.bam,
                sample_id = sample_id                               
        }     
        call MergeLogSplit {
            input:
                log_files = Trim.trimmed_log,
                sample_id = sample_id 
        }    
        call SortBam {
            input: 
                bam_file = MergeSplit.merged_bam,
                sample_id = sample_id                  
        }      
        call ByProductMetrics {
            input:
                trim_log = MergeLogSplit.merged_log,
                highconf_bam = SortBam.sorted_bam,
                sample_id = sample_id                 
        }
        call ReplaceRawReadGroup {
            input: 
                raw_bam = MergeSplit.merged_bam,
                sample_id = sample_id
        }
        call CollectInsertSizeMetrics {
            input:
                input_bam = ReplaceRawReadGroup.bam,
                sample_id = sample_id
        }
        call GroupReadByUMI {
            input:
                input_bam = ReplaceRawReadGroup.bam,
                sample_id = sample_id
        }
        if (is_captured_data) {
        call DuplexRecoveryMetrics {
            input:
                groupbyumi_bam = GroupReadByUMI.groupbyumi_bam,
                sample_id = sample_id,
                duplex_eval_bed = bait_intervals
        }
        call PlotDuplexRecoveryByTarget {
            input:
            sample_id = sample_id,
            duplex_recovery_metrics = DuplexRecoveryMetrics.duplex_recovery_metrics
        }
        }
        call FgbioCollapseReadFamilies {
            input:
                grouped_umi_bam = GroupReadByUMI.groupbyumi_bam,
                sample_id = sample_id
        }
        call AlignMolecularConsensusReads {
            input:
                mol_consensus_bam = FgbioCollapseReadFamilies.mol_consensus_bam,
                sample_id = sample_id,
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                reference_dict = reference_dict,
                reference_pac = reference_pac,
                reference_amb = reference_amb,
                reference_ann = reference_ann,
                reference_bwt = reference_bwt,
                reference_sa = reference_sa
        }
        call MergeAndSortMoleculeConsensusReads {
            input:
                mapped_sam = AlignMolecularConsensusReads.aligned_bam,
                unmapped_bam = FgbioCollapseReadFamilies.mol_consensus_bam,
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                reference_dict = reference_dict,
                sample_id = sample_id,
                sort_memory = sort_memory
        }
        if (is_captured_data) {
        call CollectSelectionMetrics {
            input: 
                reference = reference_fasta,
                reference_index = reference_fasta_index,
                reference_dict = reference_dict,
                bam_file = MergeAndSortMoleculeConsensusReads.bam,
                bam_index = MergeAndSortMoleculeConsensusReads.bai,
                sample_id = sample_id,
                bait_intervals = bait_intervals,
                target_intervals = target_intervals
        }
        call Hs_metrics {
            input:
            output_selection_metrics = CollectSelectionMetrics.output_selection_metrics
        }
        }
        call CSS_SFC_ErrorMetrics {
            input:
                ConsensusAlignedBam = MergeAndSortMoleculeConsensusReads.bam,
                ConsensusAlignedBai = MergeAndSortMoleculeConsensusReads.bai,
                sample_id = sample_id,
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                reference_dict = reference_dict,
                reference_pac = reference_pac,
                reference_amb = reference_amb,
                reference_ann = reference_ann,
                reference_bwt = reference_bwt,
                reference_sa = reference_sa,
                germline_bam = germline_bam,
                germline_bam_index = germline_bam_index,
                eval_genome_bed = eval_genome_bed
        }
        if (create_vcf) {
        call codec2MAF {
            input:
                variants_called = CSS_SFC_ErrorMetrics.variants_called,
                sample_id = sample_id
        }
        call Maf2Vcf {
            input:
            maf = codec2MAF.maf,
            ref_fasta = reference_fasta,
            ref_fai = reference_fasta_index
        }
        }
        call QC_metrics {
            input:
                byproduct_metrics = ByProductMetrics.byproduct_metrics,
                umiHistogram = GroupReadByUMI.umi_histogram,
                InsertSizeMetrics = CollectInsertSizeMetrics.insert_size_metrics,
                mutant_metrics = CSS_SFC_ErrorMetrics.mutant_metrics,
        }
        call EvalGenomeBases {
             input: 
                eval_genome_interval = eval_genome_interval
        }
        call CalculateDuplexDepth {
             input: 
              eval_genome_bases = EvalGenomeBases.eval_genome_bases,
              n_bases_eval = QC_metrics.n_bases_eval,
              n_total_fastq = QC_metrics.n_total_fastq
        }

    output {
        File byproduct_metrics = ByProductMetrics.byproduct_metrics
        File RAW_BAM = MergeSplit.merged_bam
        File RAW_BAM_index = MergeSplit.merged_bai
        File MolConsensusBAM = MergeAndSortMoleculeConsensusReads.bam
        File MolConsensusBAM_index = MergeAndSortMoleculeConsensusReads.bai
        File InsertSizeMetrics = CollectInsertSizeMetrics.insert_size_metrics
        File InsertSizeHistogram = CollectInsertSizeMetrics.insert_size_histogram
        File umiHistogram = GroupReadByUMI.umi_histogram

        File? output_selection_metrics = CollectSelectionMetrics.output_selection_metrics
        File? output_per_target_selection_metrics = CollectSelectionMetrics.output_per_target_selection_metrics
        File? output_theoretical_sensitivity = CollectSelectionMetrics.output_theoretical_sensitivity

        File mutant_metrics = CSS_SFC_ErrorMetrics.mutant_metrics
        File context_count = CSS_SFC_ErrorMetrics.context_count
        File variants_called = CSS_SFC_ErrorMetrics.variants_called

        Int n_total_fastq = QC_metrics.n_total_fastq 
        Int n_correct_products = QC_metrics.n_correct_products
        Float pct_correct_products = QC_metrics.pct_correct_products
        Int n_double_ligation = QC_metrics.n_double_ligation
        Float pct_double_ligation = QC_metrics.pct_double_ligation
        Int n_intermol = QC_metrics.n_intermol
        Float pct_intermol = QC_metrics.pct_intermol
        Int n_adp_dimer = QC_metrics.n_adp_dimer
        Float pct_adp_dimer = QC_metrics.pct_adp_dimer
        Float duplication_rate = QC_metrics.duplication_rate
        Float mean_insert_size = QC_metrics.mean_insert_size
        Int median_insert_size = QC_metrics.median_insert_size
        Int n_snv = QC_metrics.n_snv
        Int n_indel = QC_metrics.n_indel
        String n_bases_eval = QC_metrics.n_bases_eval   
        Float snv_rate = QC_metrics.snv_rate
        Float indel_rate = QC_metrics.indel_rate
        String eval_genome_bases = EvalGenomeBases.eval_genome_bases
        Float? pct_selected_bases = Hs_metrics.pct_selected_bases
        Float? mean_bait_coverage = Hs_metrics.mean_bait_coverage
        Float? mean_target_coverage = Hs_metrics.mean_target_coverage
        Float duplex_depth = CalculateDuplexDepth.duplex_depth
        Float duplex_efficiency = CalculateDuplexDepth.duplex_efficiency
        File? duplex_recovery_metrics = DuplexRecoveryMetrics.duplex_recovery_metrics
        File? ds_duplex_plots = PlotDuplexRecoveryByTarget.ds_duplex_plots
        File? outputMAF = codec2MAF.output_maf
        File? outputVCF = Maf2Vcf.vcf
    }
}

task SplitFastq1 {
    input {
        File fastq_read1
        String sample_id
        Int nsplit
        Int memory = 64
        Int? extra_disk
        Int disk_size = ceil(size(fastq_read1, "GB") * 15) + select_first([extra_disk, 0])
        String? docker_override
    }

    command <<<
        set -e
        
        zcat ~{fastq_read1} | /CODECsuite/snakemake/script/fastqsplit.pl ~{sample_id}_split_r1 ~{nsplit}

    >>>

    output {
        Array[File] split_read1 = glob("~{sample_id}_split_r1.*.fastq")
    }

    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"])
        memory: memory + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task SplitFastq2 {
    input {
        File fastq_read2
        Int nsplit
        String sample_id
        Int memory = 64
        Int? extra_disk
        Int disk_size = ceil(size(fastq_read2, "GB") * 15) + select_first([extra_disk, 0])
        String? docker_override
    }

    command <<<
        set -e
        zcat ~{fastq_read2} | /CODECsuite/snakemake/script/fastqsplit.pl ~{sample_id}_split_r2 ~{nsplit} 
        
    >>>

    output {
        Array[File] split_read2 = glob("~{sample_id}_split_r2.*.fastq")
    }

    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"])
        memory: memory + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task Trim {
    input {
        File read1
        File read2
        String sample_id
        String output_prefix
        Int split
        Int mem = 16
        Int? extra_disk
        Int disk_size = ceil(size(read1, "GB") * 8) + select_first([extra_disk, 0])
        String? docker_override
    }
        
    command {
        set -e 

        /CODECsuite/build/codec trim -1 ~{read1} -2 ~{read2} -o ~{output_prefix} -u 3 -U 3 -f 2 -t 2 -s ~{sample_id} > ~{output_prefix}.trim.log
    
    }
    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"]) 
        disks: "local-disk " + disk_size + " HDD"
        memory: mem + " GB"
    }
    output {
        File trimmed_bam = "${sample_id}_split.${split}.trim.bam"
        File trimmed_log = "${sample_id}_split.${split}.trim.log"
    }
}

task AlignRawTrimmed {
    input {
        File bam_input
        File reference_fasta
        File reference_fasta_index
        File reference_dict
        File reference_pac
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_sa
        Int mem = 16
        Int? extra_disk
        Int disk_size = ceil(size(bam_input, "GB") * 20) + select_first([extra_disk, 0])
        String sample_id
        Int split
        String? docker_override
    }

    String output_bam_name = "${sample_id}_split.${split}.aligned_tmp.bam"

    command {
        samtools fastq ~{bam_input} | \
        bwa mem \
            -K 100000000 \
            -p \
            -Y \
            ~{reference_fasta} - | samtools view - -S -b -o ~{output_bam_name} 
    }

    output {
        File aligned_bam = output_bam_name
    }

    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"])
        memory: mem + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
    }
}

task ZipperBamAlignment {
    input {
        File mapped_bam
        File unmapped_bam
        File reference_fasta
        File reference_fasta_index
        File reference_dict
        Int? mem
        Int? disk_size
        String sample_id
        Int split
        String sort_memory
        String? docker_override

    }

    String bamOutput = "${sample_id}_split.${split}.trim.aligned.bam"
    String baiOutput = "${sample_id}_split.${split}.trim.aligned.bam.bai"

    command {
        java -jar /dependencies/fgbio-2.0.2.jar --compression 0 --async-io ZipperBams \
            -i ~{mapped_bam} \
            --unmapped ~{unmapped_bam} \
            --ref ~{reference_fasta} \
        | samtools sort -m ~{sort_memory} - -o ~{bamOutput} -O BAM && samtools index ~{bamOutput}
    }

    output {
        File bam = bamOutput
        File bai = baiOutput
    }

    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"]) 
        memory: select_first([mem, 8]) + " GB"
        disks: "local-disk " + select_first([disk_size, 16]) + " HDD"
        preemptible: 3
    }
}

task MergeSplit {
    input {
        Array[File] bam_files
        String sample_id
        Int memory = 64
        Int disk_size =200
        String? docker_override
    }

    command {
        set -e
        samtools merge -@ 4 ~{sample_id}.raw.aligned.bam ~{sep=' ' bam_files} && \
        samtools index ~{sample_id}.raw.aligned.bam
    }

    output {
        File merged_bam = "~{sample_id}.raw.aligned.bam"
        File merged_bai = "~{sample_id}.raw.aligned.bam.bai"
    }

    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"]) 
        disks: "local-disk " + disk_size + " HDD"
        memory: memory + " GB"
    }
}

task MergeLogSplit {
    input {
        Array[File] log_files
        String sample_id
        Int mem = 32
        Int disk_size = 64
        String? docker_override
    }

    command {
        set -e
        python3 /CODECsuite/snakemake/script/agg_log.py ~{sep=' ' log_files} ~{sample_id}.trim.log
    }

    output {
        File merged_log = "~{sample_id}.trim.log"
    }

    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"]) 
        disks: "local-disk " + disk_size + " HDD"
        memory: mem + " GB"
    }
}

task CollectSelectionMetrics {
   input{
        File reference
        File reference_index
        File reference_dict
        File bam_file
        File bam_index
        String sample_id
        File? bait_intervals
        File? target_intervals
        Int? preemptible_attempts
        Int memory = 32
        Int disk_pad = 0
        Int disk_size = ceil(size(bam_file, "GB") * 2) + disk_pad
        String? docker_override
    }

   command {

      java -jar /dependencies/picard.jar CollectHsMetrics \
        I=~{bam_file} \
        O=~{sample_id}.selection_metrics \
        PER_TARGET_COVERAGE=~{sample_id}.per_target_selection_metrics \
        THEORETICAL_SENSITIVITY_OUTPUT=~{sample_id}.theoretical_sensitivity \
        BAIT_INTERVALS=~{bait_intervals} \
        TARGET_INTERVALS=~{target_intervals} \
        REFERENCE_SEQUENCE=~{reference}
    }

    output {
      File output_selection_metrics = "~{sample_id}.selection_metrics"
      File output_per_target_selection_metrics = "~{sample_id}.per_target_selection_metrics"
      File output_theoretical_sensitivity = "~{sample_id}.theoretical_sensitivity"
   }

   runtime {
      docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"]) 
      disks: "local-disk " + disk_size + " HDD"
      memory: memory + " GB"
      preemptible: select_first([preemptible_attempts, 3])
   }
}

task Hs_metrics {
    input {
        File output_selection_metrics
    }

    command <<<
        python3 <<CODE

        with open("~{output_selection_metrics}", 'r') as file:
            lines = file.readlines()

        metrics_line = None
        for i, line in enumerate(lines):
            if '## METRICS CLASS' in line:
                metrics_line = lines[i + 2].strip()
                break

        # Extract the values for PCT_SELECTED_BASES, MEAN_BAIT_COVERAGE, and MEAN_TARGET_COVERAGE
        if metrics_line:
            values = metrics_line.split()
            pct_selected_bases = values[6] 
            mean_bait_coverage = values[9]
            mean_target_coverage = values[32]

            with open('pct_selected_bases.txt', 'w') as f:
                f.write(pct_selected_bases + '\n')

            with open('mean_bait_coverage.txt', 'w') as f:
                f.write(mean_bait_coverage + '\n')

            with open('mean_target_coverage.txt', 'w') as f:
                f.write(mean_target_coverage + '\n')

        CODE
    >>>

    output {
      Float pct_selected_bases = read_float("pct_selected_bases.txt")     
      Float mean_bait_coverage = read_float("mean_bait_coverage.txt")
      Float mean_target_coverage = read_float("mean_target_coverage.txt")
    }

    runtime {
        memory: "16 GB"
        docker: "us.gcr.io/tag-public/metadata_upload" 
        disks: "local-disk 16 HDD"
        preemptible: 1
    }
}

task SortBam {
    input {
        File bam_file
        String sample_id
        Int mem = 64
        Int disk_size = 200
        String? docker_override
        Int? preemptible_attempts
    }

    command {
        samtools sort -n ~{bam_file} -o ~{sample_id}.raw.aligned.sortbyname.bam
    }

    output {
        File sorted_bam = "~{sample_id}.raw.aligned.sortbyname.bam"
    }

    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"]) 
        disks: "local-disk " + disk_size + " HDD"
        memory: mem + " GB"
        preemptible: select_first([preemptible_attempts, 0])
    }
}

task ByProductMetrics {
    input {
        File trim_log
        File highconf_bam
        String sample_id
        Int mem = 32
        Int disk_size = 100
        String? docker_override
    }

    command {
        python3 /CODECsuite/snakemake/script/cds_summarize.py --sample_id ~{sample_id} --trim_log ~{trim_log} \
        --highconf_bam ~{highconf_bam} > ~{sample_id}.byproduct.txt
    }

    output {
        File byproduct_metrics = "~{sample_id}.byproduct.txt"
    }

    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"]) 
        disks: "local-disk " + disk_size + " HDD"
        memory: mem + " GB"
    }
}

task ReplaceRawReadGroup {
    input {
        File raw_bam
        String sample_id
        Int memory = 64
        Int disk_size = 200
        String? docker_override
    }

    command {
        java -jar /dependencies/picard.jar AddOrReplaceReadGroups \
            I=~{raw_bam} \
            O=~{sample_id}.raw.replacerg.bam \
            CREATE_INDEX=true \
            RGID=4 \
            RGLB=lib1 \
            RGPL=ILLUMINA \
            RGPU=unit1 \
            RGSM=~{sample_id}
    }

    output {
        File bam = "~{sample_id}.raw.replacerg.bam"
        File bai = "~{sample_id}.raw.replacerg.bai"
    }

    runtime {
        memory: memory + " GB"
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"]) 
        disks: "local-disk " + disk_size + " HDD"
    }
}

task CollectInsertSizeMetrics {
    input {
        File input_bam
        String sample_id
        Int memory = 32
        Int disk_size = 200
        String? docker_override
    }

    command {
        java -jar /dependencies/picard.jar CollectInsertSizeMetrics \
            I=~{input_bam} \
            O=~{sample_id}.raw.insert_size_metrics.txt \
            H=~{sample_id}.raw.insert_size_histogram.pdf \
            M=0.5 W=600 DEVIATIONS=100
    }

    output {
        File insert_size_metrics = "~{sample_id}.raw.insert_size_metrics.txt"
        File insert_size_histogram = "~{sample_id}.raw.insert_size_histogram.pdf"
    }

    runtime {
        memory: memory + " GB"
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"])
        disks: "local-disk " + disk_size + " HDD"
    }
}

task GroupReadByUMI {
    input {
        File input_bam
        String sample_id
        Int memory = 64
        Int? extra_disk
        Int disk_size = ceil(size(input_bam, "GB") * 15) + select_first([extra_disk, 0])
        String? docker_override
        Int? preemptible_attempts
    }

    command {
        java -jar /dependencies/fgbio-2.0.2.jar --compression 1 --async-io \
            GroupReadsByUmi \
            -i ~{input_bam} \
            -o ~{sample_id}.GroupedByUmi.bam \
            -f ~{sample_id}.umiHistogram.txt \
            -m 0 \
            --strategy=paired
    }

    output {
        File groupbyumi_bam = "~{sample_id}.GroupedByUmi.bam"
        File umi_histogram = "~{sample_id}.umiHistogram.txt"
    }

    runtime {
        memory: memory + " GB"
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"])
        disks: "local-disk " + disk_size + " HDD"
        preemptible: select_first([preemptible_attempts, 0])
    }
}


task DuplexRecoveryMetrics {
    input {
        File groupbyumi_bam 
        String sample_id
        File? duplex_eval_bed
        Int memory = 32
        Int extra_disk = 0
        Int disk_size = ceil(size(groupbyumi_bam , "GB") * 3) + select_first([extra_disk, 0])
        String sorted_bam = basename(groupbyumi_bam, ".bam") + ".sorted.bam"
        String? docker_override
    }

    command {
        
        samtools sort ~{groupbyumi_bam} -o ~{sorted_bam}

        python3 /scripts/collect_duplex_metrics.py \
            --bam_file ~{sorted_bam} \
            --output_file ~{sample_id}.duplex_metrics.txt \
            --interval_list ~{duplex_eval_bed} \
            --per_target \
            --is_cds
    
    }

    output {
        File duplex_recovery_metrics = "~{sample_id}.duplex_metrics.txt"
    }

    runtime {
        memory: memory + " GB"
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"])
        disks: "local-disk " + disk_size + " HDD"
    }
}

task PlotDuplexRecoveryByTarget {
    input {
        File? duplex_recovery_metrics
        String sample_id
    }

    command <<<
        python3 <<CODE

        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns
        from matplotlib.backends.backend_pdf import PdfPages

        file_path = "~{duplex_recovery_metrics}"
        duplex_metrics = pd.read_csv(file_path, sep="\t")


        with PdfPages("~{sample_id}.ds_duplex_plots.pdf") as pdf:
            # First plot: DS Duplexes vs Read Pairs by Target
            plt.figure(figsize=(10, 6))
            for target, group in duplex_metrics.groupby('target'):
                plt.plot(group['read_pairs'], group['ds_duplexes'], label=target, linewidth = 1)
            plt.xlabel('Read Pairs')
            plt.ylabel('DS Duplexes')
            plt.title('DS Duplexes vs Read Pairs by Target')
            # plt.legend(title='Target', bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            pdf.savefig()
            plt.close()

            # Second plot: Distribution of DS Duplexes for Fraction == 1
            filtered_duplex_metrics = duplex_metrics[duplex_metrics['fraction'] == 1]
            plt.figure(figsize=(6, 6))
            sns.boxplot(data=filtered_duplex_metrics, y='ds_duplexes')
            plt.xticks([0], ['All Targets'])  # Set a single X-tick labeled "All Targets"
            plt.xlabel('')
            plt.ylabel('DS Duplexes')
            plt.title('Distribution of DS Duplexes(All Targets)')
            plt.tight_layout()
            pdf.savefig() 
            plt.close()

            print("Plots have been saved to '~{sample_id}.ds_duplex_plots.pdf'")

        CODE
    >>>


    output {
        File ds_duplex_plots = "~{sample_id}.ds_duplex_plots.pdf"
    }

    runtime {
        memory: "16 GB"
        docker: "us.gcr.io/tag-public/python:test"
        disks: "local-disk 16 HDD"
    }
}

task FgbioCollapseReadFamilies {
    input {
        File grouped_umi_bam
        String sample_id
        Int memory = 64
        Int? extra_disk
        Int disk_size = ceil(size(grouped_umi_bam, "GB") * 10) + select_first([extra_disk, 0])
        String? docker_override
        Int? preemptible_attempts        
    }

    command {
        java -jar /dependencies/fgbio-2.0.2.jar --compression 1 CallMolecularConsensusReads \
            -i ~{grouped_umi_bam} \
            -o ~{sample_id}.mol_consensus.bam \
            -p ~{sample_id} \
            --threads 2 \
            --consensus-call-overlapping-bases false \
            -M 1
    }

    output {
        File mol_consensus_bam = "~{sample_id}.mol_consensus.bam"
    }

    runtime {
        memory: memory + " GB"
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"])
        disks: "local-disk " + disk_size + " HDD"
        preemptible: select_first([preemptible_attempts, 0])
    }
}

task AlignMolecularConsensusReads {
    input {
        File mol_consensus_bam
        String sample_id
        File reference_fasta
        File reference_fasta_index
        File reference_dict
        File reference_pac
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_sa
        Int memory = 64
        Int disk_size = 200
        Int threads = 4
        Int cpu_cores = 1
        String? docker_override
        Int? preemptible_attempts
    }
        String output_bam_name = "${sample_id}.mol_consensus.aligned_tmp.bam"

    command {
        samtools fastq ~{mol_consensus_bam} \
        | bwa mem -K 100000000 -t ~{threads} -p -Y ~{reference_fasta} - |  samtools view - -S -b -o ~{output_bam_name} 
    }

    output {
        File aligned_bam = "~{sample_id}.mol_consensus.aligned_tmp.bam"
    }

    runtime {
        memory: memory + " GB"
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"]) 
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu_cores
        preemptible: select_first([preemptible_attempts, 0])
    }
}

task MergeAndSortMoleculeConsensusReads {
    input {
        File mapped_sam
        File unmapped_bam
        File reference_fasta
        File reference_fasta_index
        File reference_dict
        String sample_id
        Int memory = 64
        Int disk_size = 200
        String sort_memory
        String? docker_override
        Int? preemptible_attempts
    }

    command {
        java -jar /dependencies/fgbio-2.0.2.jar --compression 0 --async-io ZipperBams \
            -i ~{mapped_sam} \
            --unmapped ~{unmapped_bam} \
            --ref ~{reference_fasta} \
            --tags-to-reverse Consensus \
            --tags-to-revcomp Consensus \
        | samtools sort -m ~{sort_memory} - -o ~{sample_id}.mol_consensus.aligned.bam -O BAM -@ 4 && samtools index ~{sample_id}.mol_consensus.aligned.bam
    }

    output {
        File bam = "~{sample_id}.mol_consensus.aligned.bam"
        File bai = "~{sample_id}.mol_consensus.aligned.bam.bai"
    }

    runtime {
        memory: memory+ " GB"
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"]) 
        disks: "local-disk " + disk_size + " HDD"
        preemptible: select_first([preemptible_attempts, 0])
    }
}


task CSS_SFC_ErrorMetrics {
    input {
        File ConsensusAlignedBam
        File ConsensusAlignedBai
        String sample_id
        File reference_fasta
        File reference_fasta_index
        File reference_dict
        File reference_pac
        File reference_amb
        File reference_ann
        File reference_bwt
        File reference_sa
        File? germline_bam
        File? germline_bam_index
        File eval_genome_bed
        File population_based_vcf = "gs://gptag/CODEC/alfa_all.freq.breakmulti.hg38.af0001.vcf.gz"
        File population_based_vcf_index = "gs://gptag/CODEC/alfa_all.freq.breakmulti.hg38.af0001.vcf.gz.tbi"
        Int memory = 64
        Int disk_size = 200
        String? docker_override
        Int? preemptible_attempts
    }

    command {
        /CODECsuite/build/codec call -b ~{ConsensusAlignedBam} \
            -L ~{eval_genome_bed} \
            -r ~{reference_fasta} \
            -m 60 \
            -q 30 \
            -d 12 \
            ~{if defined(germline_bam) then "-n " + germline_bam else ""} \
            -V ~{population_based_vcf} \
            -x 6 \
            -c 4 \
            -5 \
            -g 30 \
            -G 250 \
            -Q 0.7 \
            -B 0.6 \
            -N 0.05 \
            -Y 10 \
            -W 1 \
            -a ~{sample_id}.mutant_metrics.txt \
            -e ~{sample_id}.variants_called.txt \
            -C ~{sample_id}.context_count.txt
    }

    output {
        File mutant_metrics = "~{sample_id}.mutant_metrics.txt"
        File variants_called = "~{sample_id}.variants_called.txt"
        File context_count = "~{sample_id}.context_count.txt"
    }

    runtime {
        memory: memory + " GB"
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec:v1.1.4"]) 
        disks: "local-disk " + disk_size + " HDD"
        preemptible: select_first([preemptible_attempts, 0])
    }
}


task QC_metrics {
    input {
        File byproduct_metrics
        File umiHistogram
        File InsertSizeMetrics
        File mutant_metrics
        Int memory = 32
        Int disk_size = 16

    }
    command <<<
        set -e

        # Use umi histogram to calculate duplication rate
        cat ~{umiHistogram} | awk -F'\t' 'NR > 1 {
            numerator += ($1 - 1) * $2;
            denominator += $1 * $2;
        } END {
            if (denominator > 0) print numerator / denominator;
            else print "NA";
        }' > duplication_rate.txt

        awk 'NR==1 {for (i=1; i<=NF; i++) if ($i=="n_total") col=i} NR==2 {print $col}' ~{byproduct_metrics} > n_total_fastq.txt
        awk 'NR==1 {for (i=1; i<=NF; i++) if ($i=="n_correct") col=i} NR==2 {print $col}' ~{byproduct_metrics} > n_correct.txt
        awk 'NR==1 {for (i=1; i<=NF; i++) if ($i=="pct_correct") col=i} NR==2 {print $col}' ~{byproduct_metrics} > pct_correct.txt
        awk 'NR==1 {for (i=1; i<=NF; i++) if ($i=="n_double_ligation") col=i} NR==2 {print $col}' ~{byproduct_metrics} > n_double_ligation.txt
        awk 'NR==1 {for (i=1; i<=NF; i++) if ($i=="pct_double_ligation") col=i} NR==2 {print $col}' ~{byproduct_metrics} > pct_double_ligation.txt
        awk 'NR==1 {for (i=1; i<=NF; i++) if ($i=="n_intermol") col=i} NR==2 {print $col}' ~{byproduct_metrics} > n_intermol.txt
        awk 'NR==1 {for (i=1; i<=NF; i++) if ($i=="pct_intermol") col=i} NR==2 {print $col}' ~{byproduct_metrics} > pct_intermol.txt
        awk 'NR==1 {for (i=1; i<=NF; i++) if ($i=="n_adp_dimer") col=i} NR==2 {print $col}' ~{byproduct_metrics} > n_adp_dimer.txt
        awk 'NR==1 {for (i=1; i<=NF; i++) if ($i=="pct_adp_dimer") col=i} NR==2 {print $col}' ~{byproduct_metrics} > pct_adp_dimer.txt

        cat ~{InsertSizeMetrics} | grep -v "#" | awk 'NR==3 {print $1}' > median_insert_size.txt
        cat ~{InsertSizeMetrics} | grep -v "#" | awk 'NR==3 {print $6}' > mean_insert_size.txt

        awk 'NR==1 {for (i=1; i<=NF; i++) if ($i=="n_snv") col=i} NR==2 {print $col}' ~{mutant_metrics} > n_snv.txt
        awk 'NR==1 {for (i=1; i<=NF; i++) if ($i=="n_indel") col=i} NR==2 {print $col}' ~{mutant_metrics} > n_indel.txt
        awk 'NR==1 {for (i=1; i<=NF; i++) if ($i=="n_bases_eval") col=i} NR==2 {print $col}' ~{mutant_metrics} > n_bases_eval.txt
        awk 'NR==1 {for (i=1; i<=NF; i++) if ($i=="snv_rate") col=i} NR==2 {print $col}' ~{mutant_metrics} > snv_rate.txt
        awk 'NR==1 {for (i=1; i<=NF; i++) if ($i=="indel_rate") col=i} NR==2 {print $col}' ~{mutant_metrics} > indel_rate.txt



    >>>
    output {
      Float duplication_rate = read_float("duplication_rate.txt")
      Int n_total_fastq = read_int("n_total_fastq.txt")
      Int n_correct_products = read_int("n_correct.txt")
      Float pct_correct_products = read_float("pct_correct.txt")
      Int n_double_ligation = read_int("n_double_ligation.txt")
      Float pct_double_ligation = read_float("pct_double_ligation.txt")
      Int n_intermol = read_int("n_intermol.txt")
      Float pct_intermol = read_float("pct_intermol.txt")
      Int n_adp_dimer = read_int("n_adp_dimer.txt")
      Float pct_adp_dimer = read_float("pct_adp_dimer.txt")
      Float mean_insert_size = read_float("mean_insert_size.txt")
      Int median_insert_size = read_int("median_insert_size.txt")
      Int n_snv = read_int("n_snv.txt")
      Int n_indel = read_int("n_indel.txt")
      String n_bases_eval = read_string("n_bases_eval.txt")     
      Float snv_rate = read_float("snv_rate.txt")
      Float indel_rate = read_float("indel_rate.txt")
    }
    runtime {
        memory: memory + " GB"
        docker: "us.gcr.io/tag-public/metadata_upload" 
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
    }
}


task EvalGenomeBases {
    input {
        File eval_genome_interval
        Int memory = 16
        Int disk_size = 8
    }

    command {
        java -jar /dependencies/picard.jar IntervalListTools \
        I=~{eval_genome_interval} COUNT_OUTPUT=eval_genome_bases.txt OUTPUT_VALUE=BASES
    
    }

    output {
        String eval_genome_bases = read_string("eval_genome_bases.txt")
    }
    runtime {
        memory: memory + " GB"
        docker: "us.gcr.io/tag-public/metadata_upload" 
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
    }
}


task CalculateDuplexDepth {
    input {
        String eval_genome_bases
        String n_bases_eval
        Int read_length = 161
        Int n_total_fastq
        Int memory = 32
        Int disk_size = 64
    }

    command <<<  
        python3 <<CODE
        
        eval_genome_bases = int("~{eval_genome_bases}")
        n_bases_eval = int("~{n_bases_eval}")
        n_total_fastq = int("~{n_total_fastq}")
        read_length = int("~{read_length}")

        duplex_depth = round (n_bases_eval / eval_genome_bases, 2)
        duplex_efficiency = round (n_bases_eval / (n_total_fastq * read_length * 2) , 4)

        with open('duplex_depth.txt', 'w') as f:
            f.write(str(duplex_depth))
        with open('duplex_efficiency.txt', 'w') as f:
            f.write(str(duplex_efficiency))

        CODE
    >>>

    output {
        Float duplex_depth = read_float('duplex_depth.txt')
        Float duplex_efficiency = read_float('duplex_efficiency.txt')
    }
    runtime {
        memory: memory + " GB"
        docker: "us.gcr.io/tag-public/metadata_upload" 
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
    }
}

task codec2MAF {
    input {
        File variants_called
        String sample_id
        String? docker_override
    }

    command <<<

    # Only run this script when there are variants called.

    if [[ $(wc -l < ~{variants_called}) -gt 2 ]]; then        
        Rscript /scripts/codec2maf.R --inputmaf ~{variants_called} --outputmaf ~{sample_id}.maf
    else 
        echo "Warning: File ~{variants_called} has no variants. Skipping this task and returning empty maf."
        touch ~{sample_id}.maf
    fi
  
    >>>

    output {
        File output_maf = "~{sample_id}.maf"
    }

    runtime{
        memory: "16G"
        docker: select_first([docker_override, "us.gcr.io/tag-public/codec2maf_r_docker:v1"])
        disks: "local-disk 16 HDD"
        preemptible: 1
    }
}


task Maf2Vcf {
  input {
    File maf
    File ref_fasta
    File ref_fai
    String outdir = "vcfs"
  }

  command <<<
    set -euo pipefail
    mkdir -p "~{outdir}"

    perl /app/script.pl \
      --input-maf "~{maf}" \
      --output-dir "~{outdir}" \
      --ref-fasta "~{ref_fasta}" \
      --ref-fai "~{ref_fai}"
  >>>

  output {
    File vcf = glob("~{outdir}/*.vcf")[0]
  }

  runtime {
    docker: "us.gcr.io/tag-public/perl-tools"
    cpu: 1
    memory: "8G"
    disks: "local-disk 16 HDD"
  }

  meta {
    description: "Convert a MAF to a multi-sample VCF using maf2vcf.pl"
  }

  parameter_meta {
    maf:        "Input MAF file"
    ref_fasta:  "Reference FASTA file"
    ref_fai:    "Reference FASTA index file (.fai)"
    outdir:     "Output directory (default: vcfs)"
  }
}
