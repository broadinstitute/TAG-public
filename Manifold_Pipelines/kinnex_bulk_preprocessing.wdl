version 1.0

workflow kinnex_bulk_preprocessing {
    meta {
        description: "End-to-end bulk MASseq pipeline: Skera split + QC plots -> Lima demux + Refine -> Merge replicates"
    }

    input {
        # --- Skera inputs ---
        File   input_bam
        File   mas_adapters_fasta
        Int    arraysize
        String? sample_id

        # --- Lima / Demux inputs ---
        File    bulk_barcodes_fasta
        Boolean trimPolyA
        Boolean clipAdapters

        # --- Merge inputs ---
        File    barcode_to_sample
        String? datasetId
        Boolean mergeBams

        # --- Shared ---
        String gcs_output_dir
        Int    num_threads = 16
    }

    # Step 1: Skera split + QC plots
    call pbSkerawQC {
        input:
            hifi_bam           = input_bam,
            sample_id          = sample_id,
            arraysize          = arraysize,
            mas_adapters_fasta = mas_adapters_fasta,
            num_threads        = num_threads,
            gcs_output_dir     = gcs_output_dir
    }

    # Step 2: Lima demux + Isoseq Refine
    call pbLimaBulk {
        input:
            skera_bam           = pbSkerawQC.skera_bam,
            sample_id           = sample_id,
            bulk_barcodes_fasta = bulk_barcodes_fasta,
            trimPolyA           = trimPolyA,
            clipAdapters        = clipAdapters,
            num_threads         = num_threads,
            gcs_output_dir      = gcs_output_dir
    }

    # Step 3: Merge replicates
    call bulkMerge {
        input:
            refine_bampath    = pbLimaBulk.refine_out,
            lima_dir          = pbLimaBulk.lima_out,
            barcode_to_sample = barcode_to_sample,
            datasetId         = datasetId,
            mergeBams         = mergeBams,
            num_threads       = num_threads,
            gcs_output_dir    = gcs_output_dir
    }

    output {
        # From Step 1
        File   skera_bam  = pbSkerawQC.skera_bam
        String QC_plots   = pbSkerawQC.QC_plots

        # From Step 2
        String refine_out = pbLimaBulk.refine_out
        String lima_out   = pbLimaBulk.lima_out

        # From Step 3
        String merge_out  = bulkMerge.merge_out
    }
}

task pbSkerawQC {
    meta {
        description: "Given hifi reads, spilts MAS 8x or 16x array structure using provided adapters"
    }
    # ------------------------------------------------
    #Inputs required
    input {
        # Required:
        File hifi_bam
        String? sample_id
        File mas_adapters_fasta
        Int num_threads
        Int arraysize = 8
        String gcs_output_dir

        # Optional:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int cpu = num_threads
        Int? boot_disk_size_gb
    }
    # Computing required disk size
    Float input_files_size_gb = 2.5*(size(hifi_bam, "GiB"))
    Int default_ram = 8
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)
    Int default_boot_disk_size_gb = 25

    # Mem is in units of GB
    Int machine_mem = select_first([mem_gb,default_ram])
    String outdir = sub(sub( gcs_output_dir + "/", "/+", "/"), "s3:/", "s3://")
    String skera_id = if defined(sample_id) then sample_id else sub(basename(hifi_bam,".bam"),".hifi_reads","")
    command <<<
        set -euxo pipefail

        echo "skera split initiated.."
        echo ~{skera_id}

        skera split -j ~{num_threads} ~{hifi_bam} ~{mas_adapters_fasta} ~{skera_id}.skera.bam
        echo "Skera split completed!"

        echo "Generating QC plots.."

        python /usr/local/src/masseq_data_processing/pb_plots/plot_concat_hist.py \
        --csv ~{skera_id}.skera.read_lengths.csv \
        --arraysize ~{arraysize} \
        --output ~{skera_id}.concat_hist.png

        python /usr/local/src/masseq_data_processing/pb_plots/plot_readlen_hist.py \
        --csv ~{skera_id}.skera.read_lengths.csv \
        --arraysize ~{arraysize} \
        --output ~{skera_id}.readlen_hist.png

        python /usr/local/src/masseq_data_processing/pb_plots/plot_ligation_heatmap.py \
        --csv ~{skera_id}.skera.ligations.csv \
        --arraysize ~{arraysize} \
        --output ~{skera_id}.ligations_heatmap.png

        echo "Copying output to S3 path provided..."
        aws s3 cp . ~{outdir}skera/ --recursive --exclude "*" --include "~{skera_id}.skera.*"
        echo "Copying skera files completed!"

        echo "Copying plots to S3 path QC_plots..."
        aws s3 cp . ~{outdir}QC_plots/ --recursive --exclude "*" --include "~{skera_id}*.png"
        echo "Copying completed!"
    >>>
    # ------------------------------------------------
    # Outputs:
    output {
        # Default output file name:
        File skera_bam        = "~{skera_id}.skera.bam"
        String QC_plots       = "~{outdir}QC_plots"
    }

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/tag-public/masseq_aws:latest"
        memory: machine_mem + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 0])
        cpu: cpu
    }

}

task pbLimaBulk {
    meta {
        description: "Given deconcatenated S-reads for Bulk samples, de-multiplexes using provided primers fasta and trims PolyA tails using Isoseq Refine"
    }
    # ------------------------------------------------
    #Inputs required
    input {
        # Required:
        File skera_bam
        String? sample_id
        File bulk_barcodes_fasta
        Boolean trimPolyA = true
        Boolean clipAdapters = true
        Int num_threads
        String gcs_output_dir
        #File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

        # Optional:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }
    # Computing required disk size
    Float input_files_size_gb = 2.5*(size(skera_bam, "GiB"))
    Int default_ram = 16
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)
    Int default_boot_disk_size_gb = 50

    # Mem is in units of GB
    Int machine_mem = if defined(mem_gb) then mem_gb else default_ram
    String outdir = sub(sub( gcs_output_dir + "/", "/+", "/"), "s3:/", "s3://")
    String isoseq_cmd = if trimPolyA then "isoseq refine --require-polya" else "isoseq refine"
    String lima_cmd = if clipAdapters then "lima --isoseq --log-level INFO" else "lima --isoseq --no-clip --log-level INFO"
    String skera_id = if defined(sample_id) then sample_id else basename(skera_bam,".skera.bam")
    command <<<
        set -euxo pipefail

        echo "Running lima demux.."
        ~{lima_cmd} -j ~{num_threads} ~{skera_bam} ~{bulk_barcodes_fasta} ~{skera_id}.lima.bam
        echo "Demuxing completed."

        echo "Copying output to S3 path provided..."
        aws s3 cp . ~{outdir}lima/ --recursive --exclude "*" --include "~{skera_id}*lima*"
        echo "Copying lima files completed!"

        echo "Running Refine..."
        for i in `ls ./*_5p--3p.bam`;
        do
         echo `basename $i`
         a=`basename $i | awk -v FS='_5p--3p.bam' '{print $1}' | awk -v FS='.lima.' '{print $1"."$2}'`
        echo $a
         ~{isoseq_cmd} -j ~{num_threads} $i ~{bulk_barcodes_fasta} ./$a.refine.bam
        done
        echo "Refine completed."

        echo "Uploading refined bams..."
        aws s3 cp . ~{outdir}refine/ --recursive --exclude "*" --include "~{skera_id}*refine*"
        echo "Copying extracted FLNC reads completed!"

    >>>
    # ------------------------------------------------
    # Outputs:
    output {
        # Default output file name:
        String refine_out  = "~{outdir}refine"
        String lima_out    = "~{outdir}lima"
    }

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/tag-public/masseq_aws:latest"
        memory: machine_mem + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 0])
        cpu: select_first([cpu, 4])
    }
}

task bulkMerge {
    meta {
        description: "Given demuxed refined reads for bulk samples, collapse replicates if mergeBams set to true else plot counts"
    }
    # ------------------------------------------------
    #Inputs required
    input {
        # Required:
        String refine_bampath
        String lima_dir
        String? datasetId = "Replicates_merged"
        File barcode_to_sample
        Boolean mergeBams = false
        Int num_threads
        String gcs_output_dir
        #File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

        # Optional:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }
    # Computing required disk size
    #Float input_files_size_gb = 2.5*(size(skera_bam, "GiB"))
    Int default_ram = 16
    Int default_disk_space_gb = 500
    #Int default_disk_space_gb = ceil((2.5 * 1024 * 2) + 1024)
    Int default_boot_disk_size_gb = 50

    # Mem is in units of GB
    Int machine_mem = if defined(mem_gb) then mem_gb else default_ram
    String outdir = sub(sub( gcs_output_dir + "/", "/+", "/"), "s3:/", "s3://")
    String refinedir = sub(sub( refine_bampath + "/", "/+", "/"), "s3:/", "s3://")
    String limadir = sub(sub( lima_dir + "/", "/+", "/"), "s3:/", "s3://")

    command <<<
        set -euxo pipefail

        echo "Fetching refined bams to combine replicates..."
        aws s3 cp ~{refinedir} . --recursive --exclude "*" --include "*refine.bam"
        echo "Fetching lima counts files.."
        aws s3 cp ~{limadir} . --recursive --exclude "*" --include "*lima.counts"

        echo "plot counts and merge"

        python /usr/local/src/masseq_data_processing/pb_plots/mergeBams.py \
            -idmap ~{barcode_to_sample} \
            -bampath . \
            -limacountsdir . \
            -outdir . \
            -mergeReplicates \
            -setTitleSamplePlot ~{datasetId}

        echo "Uploading merged bams to merge dir..."
        ls -lhrt
        aws s3 cp ./merge/ ~{outdir}merge/ --recursive
        aws s3 cp readcounts_by_sample.png ~{outdir}merge/
        aws s3 cp aggregated_lima_counts_by_sample.tsv ~{outdir}merge/
        aws s3 cp lima_counts_by_moviename.tsv ~{outdir}merge/

        echo "Completed copying merged outs and QC plot. All done"

    >>>
    # ------------------------------------------------
    # Outputs:
    output {
        # Default output file name:
        String merge_out        = "~{outdir}merge"
    }

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/tag-public/masseq_aws:latest"
        memory: machine_mem + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 0])
        cpu: select_first([cpu, 4])
    }
}

task pbSingleCell {
    meta {
        description: "Given s-reads, performs all operations from lima and isoseq pipelines."
    }
    # ------------------------------------------------
    #Inputs required
    input {
        # Required:
        File skera_bam
        String sample_id
        File primer_fasta
        File barcodes_list
        String read_design
        Boolean trimPolyA = true
        Boolean clipAdapters = true
        Int num_threads
        String gcs_output_dir
        #File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

        # Optional:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int cpu = num_threads
        Int? boot_disk_size_gb
    }

     # Computing required disk size
    Float input_files_size_gb = 2.5*(size(skera_bam, "GiB"))
    Int default_ram = 32
    Int default_disk_space_gb = ceil((input_files_size_gb * 5) + 1024)
    Int default_boot_disk_size_gb = 50

    # Mem is in units of GB
    Int machine_mem = if defined(mem_gb) then mem_gb else default_ram
    String outdir = sub(sub( gcs_output_dir + "/", "/+", "/"), "s3:/", "s3://")
    String isoseq_cmd = if trimPolyA then "isoseq refine --require-polya" else "isoseq refine"
    String lima_cmd = if clipAdapters then "lima --isoseq --log-level INFO" else "lima --isoseq --no-clip --log-level INFO"
    command <<<
        set -euxo pipefail

        cp ~{skera_bam} .
        echo "Running lima demux.."
        ~{lima_cmd} -j ~{num_threads} *.skera.bam ~{primer_fasta} ~{sample_id}.lima.bam
        echo "Demuxing completed."

        echo "Copying output to S3 path provided..."
        aws s3 cp . ~{outdir}lima/ --recursive --exclude "*" --include "~{sample_id}*lima*"
        echo "Copying lima files completed!"

        echo "Running Tag, Refine and Correct..."
        for i in `ls ./*.5p--3p.bam`;
        do
           echo `basename $i`
           a=`basename $i | awk -v FS='.5p--3p.bam' '{print $1}' | awk -v FS='.lima.' '{print $1}'`
           echo $a
           echo "tagging.."
           isoseq tag --design ~{read_design} -j ~{num_threads} $i ./$a.tagged.bam
           echo "refining.."
           ~{isoseq_cmd} -j ~{num_threads} ./$a.tagged.bam ~{primer_fasta} ./$a.refine.bam
           echo "correcting.."
           isoseq correct --barcodes ~{barcodes_list} -j ~{num_threads} ./$a.refine.bam ./~{sample_id}.corrected.bam
           echo "Done for this id!"
        done
        echo "All completed."

        echo "Uploading tagged bams..."
        aws s3 cp . ~{outdir}tag/ --recursive --exclude "*" --include "*.tagged*"
        echo "Copying extracted tagged reads completed!"

        echo "Uploading refined bams..."
        aws s3 cp . ~{outdir}refine/ --recursive --exclude "*" --include "*.refine*"
        echo "Copying refined reads completed!"

        echo "Uploading corrected bams..."
        aws s3 cp . ~{outdir}correct/ --recursive --exclude "*" --include "*.corrected*"
        echo "Copying corrected reads completed!"

    >>>
    # ------------------------------------------------
    # Outputs:
    output {
        # Default output file name:
        File corrected_reads  =  "~{sample_id}.corrected.bam"

    }

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/tag-public/masseq_aws:latest"
        memory: machine_mem + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " SSD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 0])
        cpu: cpu
    }

}

task pbGroupdedup {
    meta {
        description: "Given merged, sorted bams, performs isoseq groupdededup pipelines."
    }
    # ------------------------------------------------
    #Inputs required
    input {
        # Required:
        File input_bam
        String sample_id
        Boolean keep_non_real_cells = true
        Int num_threads
        String gcs_output_dir
        #File monitoringScript = "gs://broad-dsde-methods-tbrookin/cromwell_monitoring_script2.sh"

        # Optional:
        Int? mem_gb
        Int? preemptible_attempts
        Int? max_retries
        Int? disk_space_gb
        Int cpu = num_threads
        Int? boot_disk_size_gb
    }

     # Computing required disk size
    Float input_files_size_gb = 2.5*(size(input_bam, "GiB"))
    Int default_ram = 32
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)
    Int default_boot_disk_size_gb = 50

    # Mem is in units of GB
    Int machine_mem = select_first([mem_gb,default_ram])
    String outdir = sub(sub( gcs_output_dir + "/", "/+", "/"), "s3:/", "s3://")
    String isoseq_cmd = if keep_non_real_cells then "isoseq groupdedup --keep-non-real-cells" else "isoseq groupdedup"

    command <<<
        set -euxo pipefail

        echo "Running groupdedup.."
        ~{isoseq_cmd} -j ~{num_threads} ~{input_bam} ~{sample_id}.dedup.bam
        echo "Deduping completed."

        echo "Uploading deduped bams..."
        aws s3 cp . ~{outdir}groupdedup/ --recursive --exclude "*" --include "*.dedup*"
        echo "Copying extracted tagged reads completed!"


    >>>
    # ------------------------------------------------
    # Outputs:
    output {
        # Default output file name:
        String dedup_out = "~{outdir}groupdedup"
        File deduped_bam  =  "~{sample_id}.dedup.bam"
    }

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/tag-public/masseq_aws:latest"
        memory: machine_mem + " GiB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " SSD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: select_first([preemptible_attempts, 0])
        maxRetries: select_first([max_retries,0])
        cpu: cpu
    }

}
