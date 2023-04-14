#Author: Brian Granger, Micah Rickles-Young
#Date: 4/3/20
#Snapshot 42
#This method is for taking SmartSeq2 qc output and running an R script on the results to try to provide qc at the plate level.

workflow SC_plate{
	File RPlateQC
	File metadata
	File annot_gtf
    String flowcells
    String? LCSET
    String species_name
    Array[String] smid
    Array[String]? cell_types
    Array[File] aln_list
	Array[File] base_list
    Array[File] dup_list
    Array[File] insert_list
    Array[File] rna_list
    Array[File] qual_list
    Array[File] rsem_list
	Array[File]? adapt_list
    Array[String] names

        call graphPlate {
         input:
           RPlateQC = RPlateQC,
	       metadata = metadata,
	       annot_gtf = annot_gtf,
           flowcells = flowcells,
           LCSET = LCSET,
           smid=smid,
           species_name=species_name,
           cell_types = cell_types,
           aln_list = aln_list,
	       base_list = base_list,
	       dup_list = dup_list,
	       insert_list = insert_list,
	       rna_list = rna_list,
	       qual_list = qual_list,
	       rsem_list = rsem_list,
           adapt_list = adapt_list,
           names = names
        }

		call gsutil_cp{
			input:
			plate_qc_metrics = graphPlate.plate_qc_metrics
		}
}

task graphPlate{
	File RPlateQC
	File metadata
    String metadata_basename = basename(metadata,".metadata.txt")
	File annot_gtf
    String flowcells
    String? LCSET
    String species_name
    Array[String]? cell_types
    Array[String] smid
    Array[File] aln_list
	Array[File] base_list
	Array[File] dup_list
	Array[File] insert_list
	Array[File] rna_list
	Array[File] qual_list
	Array[File] rsem_list
    Array[File]? adapt_list
    Array[String] names

    Float? extra_mem
    Float memory = 7.5 + select_first([extra_mem,0])

    Int? extra_space
    Int disk_space = 500 + select_first([extra_space,0])

    Int? extra_boot_space
    Int boot_disk_space = 10 + select_first([extra_boot_space, 0])

        command <<<
        set -euo pipefail
		# First we have to make a bunch of folders for all the different types of files that we have from the single cell qc. We need aln_sum, base_call, dup_met, insert_met, rna_cov, qual_cyc, rsem_gene
		# So time for the first one: aln_sum
        mkdir aln_sum/

		mv ${sep=" " aln_list} aln_sum/
		# Ok, all the files should be moved (and renamed, not sure I want that, but we'll see). Let's check them out.
        echo 'aln_sum:'
        ls aln_sum/

		# Second folder: base_call
        mkdir base_call/

		mv ${sep=" " base_list} base_call/

		#Ok, again, all files moved. Let's look
		echo 'base_call:'
		ls base_call/


		# Third folder: dup_met
        mkdir dup_met/

		mv ${sep=" " dup_list} dup_met/
		#Ok, again, all files moved. Let's look
		echo 'dup_met:'
		ls dup_met/

		# Fourth folder: insert_met
		mkdir insert_met/

		mv ${sep=" " insert_list} insert_met/

		echo 'insert_met:'
		ls insert_met/


		# Fifth folder: rna_cov
		mkdir rna_cov/

		mv ${sep=" " rna_list} rna_cov/

        # This folder's special. This is where we have files that may or may not have a histogram. Current understanding is that it should be all 0s.

        #create a histogram file
		echo -e "## HISTOGRAM\tjava.lang.Integer\nnormalized_position\tAll_Reads.normalized_coverage" > histo.txt

		for i in `seq 0 100`; do
			echo -e "$i\t0.0" >> histo.txt
		done
		echo "" >> histo.txt

		# add the histogram section to any file in the rna_cov/ folder that's missing it.
		for filename in rna_cov/*; do
			read lines f <<< $(wc -l $filename)

			if [ $lines -eq '10' ]
			then
				head -n 9 $filename > temp1.txt
				cat temp1.txt histo.txt > temp2.txt
				cp temp2.txt $filename
			fi
		done

		# clean up temporary files
        if [ -f temp1.txt ]; then
			rm temp1.txt
        fi
        if [ -f temp2.txt ]; then
        	rm temp2.txt
        fi

		echo 'rna_cov:'
		ls rna_cov/
		wc -l rna_cov/*

		# Sixth folder: qual_cyc
		mkdir qual_cyc/

		mv ${sep=" " qual_list} qual_cyc/

		echo 'qual_cyc:'
		ls qual_cyc/


		# Seventh folder: rsem_gene
		mkdir rsem_gene/

		mv ${sep=" " rsem_list} rsem_gene/

		echo 'rsem_gene:'
		ls rsem_gene/

		# Eighth folder: adapt_content
		mkdir adapt_content/

		if [ "${sep=" " adapt_list}" != "" ]; then
        	mv ${sep=" " adapt_list} adapt_content/
		fi
		echo 'adapt_content:'
		ls adapt_content/

		# Ok, finished moving all the files, ready to run the script after creating an images folder for the output (it might be created by the script, but let's be sure
		mkdir images

		export R_MAX_MEM_SIZE=750	#I don't think this ends up doing anything honestly.

        CELL_TYPES="$(echo ${sep="," cell_types} | sed 's/ /_/g')"
		# R-3.4.0 location                      Rscript in  f1       f2         f3       f4          f5       f6                 f7         8          9
        /usr/tag/software/R/R-3.4.0/bin/Rscript ${RPlateQC} aln_sum/ base_call/ dup_met/ insert_met/ rna_cov/ qual_cyc/ rna_cov/ rsem_gene/ adapt_content/ ${metadata} ${annot_gtf} ${species_name} ${flowcells} ${sep="," smid} $CELL_TYPES ${LCSET}

        echo "Finished running R script\n"

        # each plot is in a separate pdf. I want to combine these into 2 relevant pdfs. We're going to use ghostscript:
        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=${metadata_basename}.sequencingqc.pdf p3.pdf p7.pdf p8.pdf p1.pdf p2.pdf p13.pdf
		gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=${metadata_basename}.transcriptqc.pdf p10.pdf p9.pdf p5.pdf p11.pdf

        tar -cz images processedQC.Rdata > ${metadata_basename}.images.tar.gz

         head -n1 ${metadata_basename}.plate_qc_metrics.txt > ${metadata_basename}.plate_qc_metrics_temp.txt
         tail -n +2 ${metadata_basename}.plate_qc_metrics.txt | sort -k1,1 -k2,2n >> ${metadata_basename}.plate_qc_metrics_temp.txt
         mv ${metadata_basename}.plate_qc_metrics_temp.txt ${metadata_basename}.plate_qc_metrics.txt

        echo "Reached end of WDL"
        >>>

        output {
		# all our output is in the images folder, plus plots and Rdata in cwd
		File images = "${metadata_basename}.images.tar.gz"
        File plate_summary_metrics = "${metadata_basename}.plate_summary_metrics.txt"
        File sequence_plots = "${metadata_basename}.sequencingqc.pdf"
        File transcript_plots = "${metadata_basename}.transcriptqc.pdf"
        File plate_qc_metrics = "${metadata_basename}.plate_qc_metrics.txt"
        }

	runtime {
	    docker: "bgranger/ss2_qc:0.1"
        memory: memory + "GB"
        cpu: "2"
		disks: "local-disk "+disk_space+" HDD"
        bootDiskSizeGb: boot_disk_space
    }
}

task gsutil_cp{
	File plate_qc_metrics
	String? target_google_bucket = "gs://fc-735a9d10-0cf6-4ae5-a203-5e5522bf5c3c/tableau_files"

	command <<<
		#Run gsutil cp and capture its exit status
		gsutil cp ${plate_qc_metrics} ${target_google_bucket}
		gsutil_exit_status=$?

		# Check if gsutil cp was successful
		if [[ $gsutil_exit_status -eq 0 ]]; then
		  echo "gsutil cp succeeded"
		else
		  echo "gsutil cp failed with exit code $gsutil_exit_status"
		fi
	>>>

	runtime{
		docker: "gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
	}

}
