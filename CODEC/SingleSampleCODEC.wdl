version 1.0

workflow SingleSampleCODEC {
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
        File germline_bam
        File germline_bam_index
        Int num_parallel
        String sort_memory
        File eval_genome_interval = "gs://gptag/CODEC/GRCh38_notinalldifficultregions.interval_list"
        File eval_genome_bed = "gs://gptag/CODEC/GRCh38_notinalldifficultregions.bed"
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
        call MarkRawDuplicates {
            input:
                input_bam = ReplaceRawReadGroup.bam,
                sample_id = sample_id
        }
        call CollectInsertSizeMetrics {
            input:
                input_bam = MarkRawDuplicates.dup_marked_bam,
                sample_id = sample_id
        }
        call GroupReadByUMI {
            input:
                input_bam = MarkRawDuplicates.dup_marked_bam,
                sample_id = sample_id
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
        call CollectWgsMetrics {
            input:
                ConsensusAlignedBam = MergeAndSortMoleculeConsensusReads.bam,
                ConsensusAlignedBai = MergeAndSortMoleculeConsensusReads.bai,
                sample_id = sample_id,
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                reference_dict = reference_dict,
                eval_genome_interval = eval_genome_interval
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
        call QC_metrics {
            input:
                byproduct_metrics = ByProductMetrics.byproduct_metrics,
                WgsMetrics = CollectWgsMetrics.WgsMetrics,
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
              mean_insert_size = QC_metrics.mean_insert_size,
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
        File DuplicationMetrics = MarkRawDuplicates.dup_metrics
        File umiHistogram = GroupReadByUMI.umi_histogram
        File WgsMetrics = CollectWgsMetrics.WgsMetrics
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
        Float raw_dedupped_mean_cov = QC_metrics.raw_dedupped_mean_cov
        Int raw_dedupped_median_cov = QC_metrics.raw_dedupped_median_cov
        Float duplication_rate = QC_metrics.duplication_rate
        Float mean_insert_size = QC_metrics.mean_insert_size
        Int median_insert_size = QC_metrics.median_insert_size
        Int n_snv = QC_metrics.n_snv
        Int n_indel = QC_metrics.n_indel
        String n_bases_eval = QC_metrics.n_bases_eval   
        Float snv_rate = QC_metrics.snv_rate
        Float indel_rate = QC_metrics.indel_rate
        String eval_genome_bases = EvalGenomeBases.eval_genome_bases
        Float duplex_depth = CalculateDuplexDepth.duplex_depth
        Float duplex_efficiency = CalculateDuplexDepth.duplex_efficiency

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
    }

    command <<<
        set -e
        
        zcat ~{fastq_read1} | /CODECsuite/snakemake/script/fastqsplit.pl ~{sample_id}_split_r1 ~{nsplit}

    >>>

    output {
        Array[File] split_read1 = glob("~{sample_id}_split_r1.*.fastq")
    }

    runtime {
        docker: "us.gcr.io/tag-team-160914/codec:v1"
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
    }

    command <<<
        set -e
        zcat ~{fastq_read2} | /CODECsuite/snakemake/script/fastqsplit.pl ~{sample_id}_split_r2 ~{nsplit} 
        
    >>>

    output {
        Array[File] split_read2 = glob("~{sample_id}_split_r2.*.fastq")
    }

    runtime {
        docker: "us.gcr.io/tag-team-160914/codec:v1"
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
    }
        
    command {
        set -e 

        /CODECsuite/build/codec trim -1 ~{read1} -2 ~{read2} -o ~{output_prefix} -u 3 -U 3 -f 2 -t 2 -s ~{sample_id} > ~{output_prefix}.trim.log
    
    }
    runtime {
        docker: "us.gcr.io/tag-team-160914/codec:v1" 
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
        docker: "us.gcr.io/tag-team-160914/codec:v1"
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
        docker: "us.gcr.io/tag-team-160914/codec:v1" 
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
        docker: "us.gcr.io/tag-team-160914/codec:v1" 
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
    }

    command {
        set -e
        python3 /CODECsuite/snakemake/script/agg_log.py ~{sep=' ' log_files} ~{sample_id}.trim.log
    }

    output {
        File merged_log = "~{sample_id}.trim.log"
    }

    runtime {
        docker: "us.gcr.io/tag-team-160914/codec:v1" 
        disks: "local-disk " + disk_size + " HDD"
        memory: mem + " GB"
    }
}

task SortBam {
    input {
        File bam_file
        String sample_id
        Int mem = 64
        Int disk_size = 200
    }

    command {
        samtools sort -n ~{bam_file} -o ~{sample_id}.raw.aligned.sortbyname.bam
    }

    output {
        File sorted_bam = "~{sample_id}.raw.aligned.sortbyname.bam"
    }

    runtime {
        docker: "us.gcr.io/tag-team-160914/codec:v1" 
        disks: "local-disk " + disk_size + " HDD"
        memory: mem + " GB"
        preemptible: 2
    }
}

task ByProductMetrics {
    input {
        File trim_log
        File highconf_bam
        String sample_id
        Int mem = 32
        Int disk_size = 100
    }

    command {
        python3 /CODECsuite/snakemake/script/cds_summarize.py --sample_id ~{sample_id} --trim_log ~{trim_log} \
        --highconf_bam ~{highconf_bam} > ~{sample_id}.byproduct.txt
    }

    output {
        File byproduct_metrics = "~{sample_id}.byproduct.txt"
    }

    runtime {
        docker: "us.gcr.io/tag-team-160914/codec:v1" 
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
        docker: "us.gcr.io/tag-team-160914/codec:v1" 
        disks: "local-disk " + disk_size + " HDD"
    }
}

task MarkRawDuplicates {
    input {
        File input_bam
        String sample_id
        Int memory = 64
        Int disk_size = 200
    }

    command {
        java -jar /dependencies/picard.jar MarkDuplicates \
            I=~{input_bam} \
            O=~{sample_id}.raw.replacerg.markdup.bam \
            M=~{sample_id}.raw.marked_duplicates.txt \
            CREATE_INDEX=true \
            TAG_DUPLICATE_SET_MEMBERS=true \
            TAGGING_POLICY=All
            
        samtools index ~{sample_id}.raw.replacerg.markdup.bam
    }

    output {
        File dup_marked_bam = "~{sample_id}.raw.replacerg.markdup.bam"
        File dup_marked_bai = "~{sample_id}.raw.replacerg.markdup.bam.bai"
        File dup_metrics = "~{sample_id}.raw.marked_duplicates.txt"
    }

    runtime {
        memory: memory + " GB"
        docker: "us.gcr.io/tag-team-160914/codec:v1" 
        disks: "local-disk " + disk_size + " HDD"
    }
}

task CollectInsertSizeMetrics {
    input {
        File input_bam
        String sample_id
        Int memory = 32
        Int disk_size = 200
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
        docker: "us.gcr.io/tag-team-160914/codec:v1"
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
        docker: "us.gcr.io/tag-team-160914/codec:v1"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task FgbioCollapseReadFamilies {
    input {
        File grouped_umi_bam
        String sample_id
        Int memory = 64
        Int? extra_disk
        Int disk_size = ceil(size(grouped_umi_bam, "GB") * 10) + select_first([extra_disk, 0])
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
        docker: "us.gcr.io/tag-team-160914/codec:v1"
        disks: "local-disk " + disk_size + " HDD"
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
        docker: "us.gcr.io/tag-team-160914/codec:v1" 
        disks: "local-disk " + disk_size + " HDD"
        cpu: cpu_cores
        preemptible: 3
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
        docker: "us.gcr.io/tag-team-160914/codec:v1" 
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
    }
}

task CollectWgsMetrics {
    input {
        File ConsensusAlignedBam
        File ConsensusAlignedBai
        String sample_id
        File reference_fasta
        File reference_fasta_index
        File reference_dict
        File eval_genome_interval
        Int memory = 64
        Int disk_size = 200
    }

    command {
        java -jar /dependencies/picard.jar CollectWgsMetrics \
        I=~{ConsensusAlignedBam} O=~{sample_id}.wgs_metrics.txt R=~{reference_fasta} INTERVALS=~{eval_genome_interval} \
        COUNT_UNPAIRED=true MINIMUM_BASE_QUALITY=0 MINIMUM_MAPPING_QUALITY=0
    }

    output {
        File WgsMetrics = "~{sample_id}.wgs_metrics.txt"
    }

    runtime {
        memory: memory + " GB"
        docker: "us.gcr.io/tag-team-160914/codec:v1" 
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
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
        File germline_bam
        File germline_bam_index
        File eval_genome_bed
        File population_based_vcf = "gs://gptag/CODEC/alfa_all.freq.breakmulti.hg38.af0001.vcf.gz"
        File population_based_vcf_index = "gs://gptag/CODEC/alfa_all.freq.breakmulti.hg38.af0001.vcf.gz.tbi"
        Int memory = 64
        Int disk_size = 200
    }

    command {
        /CODECsuite/build/codec call -b ~{ConsensusAlignedBam} \
            -L ~{eval_genome_bed} \
            -r ~{reference_fasta} \
            -m 60 \
            -q 30 \
            -d 12 \
            -n ~{germline_bam} \
            -V ~{population_based_vcf} \
            -x 6 \
            -c 4 \
            -5 \
            -g 30 \
            -G 250 \
            -Q 0.7 \
            -B 0.6 \
            -N 0.05 \
            -Y 5 \
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
        docker: "us.gcr.io/tag-team-160914/codec:v1" 
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 3
    }
}


task QC_metrics {
    input {
        File byproduct_metrics
        File WgsMetrics
        File umiHistogram
        File InsertSizeMetrics
        File mutant_metrics
        Int memory = 16
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

        cat ~{WgsMetrics} | grep -v "#" | awk 'NR==3 {print $2}' > raw_dedupped_mean_cov.txt
        cat ~{WgsMetrics} | grep -v "#" | awk 'NR==3 {print $4}' > raw_dedupped_median_cov.txt
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
      Float raw_dedupped_mean_cov = read_float("raw_dedupped_mean_cov.txt")
      Int raw_dedupped_median_cov = read_int("raw_dedupped_median_cov.txt")
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
        docker: "us.gcr.io/tag-team-160914/picard_docker" 
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
        docker: "us.gcr.io/tag-team-160914/picard_docker" 
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
    }
}


task CalculateDuplexDepth {
    input {
        String eval_genome_bases
        String n_bases_eval
        Float mean_insert_size
        Int n_total_fastq
        Int memory = 16
        Int disk_size = 8
    }

    command <<<  
        python3 <<CODE
        
        eval_genome_bases = int("~{eval_genome_bases}")
        n_bases_eval = int("~{n_bases_eval}")
        mean_insert_size = float("~{mean_insert_size}")
        n_total_fastq = int("~{n_total_fastq}")

        raw_read_depth = round (n_total_fastq * mean_insert_size / eval_genome_bases, 2)
        duplex_depth = round (n_bases_eval / eval_genome_bases, 2)
        duplex_efficiency = round (duplex_depth / raw_read_depth , 4)

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
        docker: "us.gcr.io/tag-team-160914/picard_docker" 
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
    }
}