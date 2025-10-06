version 1.0

workflow DeliverCNVOutputs {
    input {
        String delivery_workspace_bucket
        String this_workspace
        String this_namespace
        String delivery_namespace
        String delivery_workspace
        Array[String] cnv_pass_fail
        Array[String] pair_name
        Array[String] called_copy_ratio_legacy_segments_normal
        Array[String] called_copy_ratio_legacy_segments_tumor
        Array[String] called_copy_ratio_segments_normal
        Array[String] called_copy_ratio_segments_tumor
        Array[String] denoised_copy_ratios_normal
        Array[String] denoised_copy_ratios_plot_normal
        Array[String] denoised_copy_ratios_plot_tumor
        Array[String] denoised_copy_ratios_tumor
        Array[String] oncotated_called_file_tumor
        Array[String] oncotated_called_gene_list_file_tumor
        Array[String] allelic_counts_tumor_files
        Array[String] allelic_counts_normal_files
        Array[String] AllChrPlot
        Array[String] QUICvizPDF
    }
    scatter (i in range(length(pair_name))){
        call CreateDeliverablesList {
            input:
            pair_name = pair_name[i],
            called_copy_ratio_legacy_segments_normal = called_copy_ratio_legacy_segments_normal[i],
            called_copy_ratio_legacy_segments_tumor = called_copy_ratio_legacy_segments_tumor[i],
            called_copy_ratio_segments_normal = called_copy_ratio_segments_normal[i],
            called_copy_ratio_segments_tumor = called_copy_ratio_segments_tumor[i],
            denoised_copy_ratios_normal = denoised_copy_ratios_normal[i],
            denoised_copy_ratios_plot_normal = denoised_copy_ratios_plot_normal[i],
            denoised_copy_ratios_plot_tumor = denoised_copy_ratios_plot_tumor[i],
            denoised_copy_ratios_tumor = denoised_copy_ratios_tumor[i],
            oncotated_called_file_tumor = oncotated_called_file_tumor[i],
            oncotated_called_gene_list_file_tumor = oncotated_called_gene_list_file_tumor[i],
            allelic_counts_tumor = allelic_counts_tumor_files[i],
            allelic_counts_normal = allelic_counts_normal_files[i],
            AllChrPlot = AllChrPlot[i],
            QUICvizPDF = QUICvizPDF[i]
        }
        if(cnv_pass_fail[i] == "PASS"){
            call DeliverCNVOutputsTask {
                input:
                pair_name = pair_name[i],
                deliverable_files_list = CreateDeliverablesList.deliverable_files_list,
                delivery_workspace_bucket = delivery_workspace_bucket
            }
        }
    }
    
    call DeliverPairTable {
        input: 
            delivery_pairs = pair_name,
            this_namespace = this_namespace,
            this_workspace = this_workspace,
            delivery_namespace = delivery_namespace,
            delivery_workspace = delivery_workspace,
            pair_delivered_messages = select_all(DeliverCNVOutputsTask.output_message)
    }

output {
        Array[Array[File]] md5sums = select_all(DeliverCNVOutputsTask.md5_files)
        String pair_table_upload_message = DeliverPairTable.pair_upload_message
    }
}

task DeliverCNVOutputsTask {
    input {
        String delivery_workspace_bucket
        String pair_name
        File deliverable_files_list
    }

    command <<<
        for file in `cat ~{deliverable_files_list}`; do
            file_name=`basename $file`
            gsutil hash -mh $file | grep md5 | cut -f4 > $file_name".md5"
            gsutil cp $file "gs://"~{delivery_workspace_bucket}"/"~{pair_name}"/"
            gsutil cp $file_name".md5" "gs://"~{delivery_workspace_bucket}"/"~{pair_name}"/"
        done
        echo "All files delivered for ~{pair_name}."
    >>>

    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        preemptible_tries:     1
        max_retries:           1
        docker:"us.gcr.io/google.com/cloudsdktool/google-cloud-cli:alpine"
    }

    output {
        Array[File] md5_files = glob("*.md5")
        String output_message = read_string(stdout())
    }
}

task CreateDeliverablesList {
    input {
        String pair_name
        String called_copy_ratio_legacy_segments_normal
        String called_copy_ratio_legacy_segments_tumor
        String called_copy_ratio_segments_normal
        String called_copy_ratio_segments_tumor
        String denoised_copy_ratios_normal
        String denoised_copy_ratios_plot_normal
        String denoised_copy_ratios_plot_tumor
        String denoised_copy_ratios_tumor
        String oncotated_called_file_tumor
        String oncotated_called_gene_list_file_tumor
        String allelic_counts_tumor
        String allelic_counts_normal
        String AllChrPlot
        String QUICvizPDF
    }

    command <<<
        touch ~{pair_name}.deliverables.txt
        echo ~{called_copy_ratio_legacy_segments_normal} >> ~{pair_name}.deliverables.txt
        echo ~{called_copy_ratio_legacy_segments_tumor} >> ~{pair_name}.deliverables.txt
        echo ~{called_copy_ratio_segments_normal} >> ~{pair_name}.deliverables.txt
        echo ~{called_copy_ratio_segments_tumor} >> ~{pair_name}.deliverables.txt
        echo ~{denoised_copy_ratios_normal} >> ~{pair_name}.deliverables.txt
        echo ~{denoised_copy_ratios_tumor} >> ~{pair_name}.deliverables.txt
        echo ~{denoised_copy_ratios_plot_normal} >> ~{pair_name}.deliverables.txt
        echo ~{denoised_copy_ratios_plot_tumor} >> ~{pair_name}.deliverables.txt
        echo ~{oncotated_called_file_tumor} >> ~{pair_name}.deliverables.txt
        echo ~{oncotated_called_gene_list_file_tumor} >> ~{pair_name}.deliverables.txt
        echo ~{allelic_counts_tumor} >> ~{pair_name}.deliverables.txt
        echo ~{allelic_counts_normal} >> ~{pair_name}.deliverables.txt
        echo ~{AllChrPlot} >> ~{pair_name}.deliverables.txt
        echo ~{QUICvizPDF} >> ~{pair_name}.deliverables.txt
    >>>

    output {
        File deliverable_files_list = "~{pair_name}.deliverables.txt"
    }

    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        preemptible_tries:     1
        max_retries:           1
        docker:"us.gcr.io/tag-team-160914/tag-tools:1.0.0"
    }
}

task DeliverPairTable {
    input {
        Array[String] delivery_pairs
        Array[String] pair_delivered_messages
        String this_namespace
        String this_workspace
        String delivery_namespace
        String delivery_workspace
    }
    command <<<
        source activate NeoVax-Input-Parser
        python3 <<CODE

        import firecloud.api as fapi
        import pandas as pd

        this_namespace = "~{this_namespace}"
        this_workspace = "~{this_workspace}"
        delivery_namespace = "~{delivery_namespace}"
        delivery_workspace = "~{delivery_workspace}"
        pairs_to_deliver = "~{sep=',' delivery_pairs}".split(",")
        delivered_pairs = pd.DataFrame(fapi.get_entities(delivery_namespace, delivery_workspace, etype='pair').json())
        if len(delivered_pairs) > 0:
            duplicate_pairs = [p for p in delivered_pairs['name'] if str(p) in pairs_to_deliver]
            assert len(duplicate_pairs) == 0, f'Attempting to deliver pairs that already exist in the delivery workspace: {",".join(duplicate_pairs)}'

        new_pairs = pd.DataFrame(fapi.get_entities(this_namespace, this_workspace, etype='pair').json())
        pair_df = pd.DataFrame(new_pairs['name'].tolist())
        pair_df.columns = ['entity:pair_id']
        expanded_attribute_df = pd.DataFrame(new_pairs['attributes'].tolist())
        new_pair_df = pd.concat([pair_df, expanded_attribute_df], axis=1)
        new_pair_df.index = new_pair_df['entity:pair_id']
        new_pair_df = new_pair_df[[p in pairs_to_deliver for p in new_pair_df.index]]
        new_pair_df = new_pair_df[['cnv_pass_fail']]
        new_pair_df.to_csv("delivery_pairs.tsv",sep="\t", index=True)
        # Upload generated tsvs to data model in analysis workspace
        upload_pair_response = fapi.upload_entities_tsv(delivery_namespace, delivery_workspace, './delivery_pairs.tsv', model='flexible')
        if "message" not in upload_pair_response:
            print("Pair table was uploaded successfully!")
        else:
            print("Error uploading pair table:", upload_pair_response["message"])
        CODE
        >>>
    runtime {
       docker: "us.gcr.io/tag-team-160914/neovax-parsley:2.2.1.0"
       preemptible: 0
    }
    output {
        String pair_upload_message = read_string(stdout())
    }
}