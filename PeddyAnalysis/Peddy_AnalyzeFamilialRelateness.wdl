version 1.0


workflow Peddy_AnalyzeFamilialRelateness {
    input {
        Array[String] sample_ids
        Array[String] family_ids
        Array[File] gvcf_paths
        Array[File] gvcf_index_paths
        Array[String] pedigrees
        Array[String] reported_sexes
        Int buffer_disk_size = 20
    }
        Int merging_disk_size = ceil(length(gvcf_paths) + buffer_disk_size)

    call AnalyzeAndCheckFamilySamples {
        input:  
        sample_ids = sample_ids, 
        family_ids = family_ids
    }
    if (AnalyzeAndCheckFamilySamples.num_single_sample_families == 0) {
        scatter (family_id in AnalyzeAndCheckFamilySamples.unique_family_ids) {
            call MergeFamilyVCFs {
                input:
                    family_id = family_id,
                    sample_ids = sample_ids,
                    gvcf_paths = gvcf_paths,
                    gvcf_index_paths = gvcf_index_paths,
                    family_ids = family_ids,
                    pedigrees = pedigrees,
                    reported_sexes = reported_sexes,
                    disk_size = merging_disk_size
            }
            call RunPlink {
                input:
                    family_id =  family_id,
                    merged_vcf = MergeFamilyVCFs.merged_vcf,
                    merged_vcf_index = MergeFamilyVCFs.merged_vcf_index
            }
            call UpdateFamFile {
                input:
                    fam_file = RunPlink.binary_fam,
                    known_trio_info = MergeFamilyVCFs.family_info,
                    family_id = family_id
            }
            call RunPeddy {
                input:
                    prefix = family_id,
                    merged_vcf = MergeFamilyVCFs.merged_vcf,
                    merged_vcf_index = MergeFamilyVCFs.merged_vcf_index,
                    fam_file = UpdateFamFile.updated_fam_file
            }
        }
        call MergePeddyResults {
            input:
                peddy_results = RunPeddy.pedigree_prediction_stats
        }

        call PlotPeddyResults {
            input:
                merged_peddy_results = MergePeddyResults.merged_peddy_results
        }
    }



    output {
        File family_composition_log =  AnalyzeAndCheckFamilySamples.family_composition_log
        Int num_single_sample_families =  AnalyzeAndCheckFamilySamples.num_single_sample_families
        Array[File]? merged_vcf_files = MergeFamilyVCFs.merged_vcf
        Array[File]? merged_vcf_indices = MergeFamilyVCFs.merged_vcf_index
        Array[File]? family_info_files = MergeFamilyVCFs.family_info
        Array[File]? updated_fam_file = UpdateFamFile.updated_fam_file
        File? merged_peddy_results = MergePeddyResults.merged_peddy_results
        File? all_family_peddy_prediction_plot = PlotPeddyResults.all_family_peddy_prediction_plot
    }
    
}


task AnalyzeAndCheckFamilySamples {
    input {
        Array[String] sample_ids
        Array[String] family_ids
        Int memory = 16
        Int disk_size = 16
    }

    command <<<
        # Write sample_ids and family_ids to file to make sure they are referred correctly in python
        echo ~{sep=' ' sample_ids} > sample_ids.txt
        echo ~{sep=' ' family_ids} > family_ids.txt

        python3 <<CODE
        import pandas as pd

        def analyze_family_samples(sample_ids, family_ids, log_file):
            df = pd.DataFrame({'entity:sample_id': sample_ids, 'sidr_family_id': family_ids})

            total_samples = len(df)
            distinct_family_ids = df['sidr_family_id'].nunique()  # Count distinct family IDs
            family_counts = df['sidr_family_id'].value_counts()  # Count the number of samples per family
            num_trios = (family_counts == 3).sum()
            num_duos = (family_counts == 2).sum()
            num_single = (family_counts == 1).sum()

            with open(log_file, 'w') as log:
                # Check for families with only one sample
                single_sample_families = family_counts[family_counts == 1]
                if not single_sample_families.empty:
                    for family_id in single_sample_families.index:
                        sample_id = df[df['sidr_family_id'] == family_id]['entity:sample_id'].values[0]
                        log.write(f'Warning: Family {family_id} has only one sample: {sample_id}\n')

                log.write(f'This sample set has {num_trios} Trios and {num_duos} Duos.\n')
                log.write(f'There are {distinct_family_ids} distinct families in total.\n')
                log.write(f'The total number of samples is {total_samples}.\n')

            return num_single, single_sample_families.index.tolist()

        if __name__ == "__main__":
            # Read sample_ids and family_ids from text files
            with open("sample_ids.txt", "r") as f:
                sample_ids = f.read().strip().split(' ')
            
            with open("family_ids.txt", "r") as f:
                family_ids = f.read().strip().split(' ')

            log_file = "analyze_family_samples.log"

            num_single, single_sample_family_ids = analyze_family_samples(sample_ids, family_ids, log_file)

            # Save unique family IDs to a file
            unique_family_ids = list(set(family_ids))
            with open('unique_family_ids.txt', 'w') as f:
                for family_id in unique_family_ids:
                    f.write(f"{family_id}\n")
            # Save the number of single sample families to a file
            with open('num_single_sample_families.txt', 'w') as f:
                f.write(f"{num_single}\n")

        CODE
        >>>

    output {
        File family_composition_log = "analyze_family_samples.log"
        Int num_single_sample_families = read_int("num_single_sample_families.txt")
        Array[String] unique_family_ids = read_lines("unique_family_ids.txt")
    }

    runtime {
        docker: "us.gcr.io/tag-public/peddy-analysis:v1"
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task MergeFamilyVCFs {
    input {
        String family_id
        Array[String] family_ids
        Array[String] sample_ids
        Array[File] gvcf_paths
        Array[File] gvcf_index_paths
        Array[String] pedigrees
        Array[String] reported_sexes
        Int memory = 32
        Int disk_size
    }

    command <<<
        # Merge and Index the multi-sample gVCF files for the family

        echo ~{sep=' ' sample_ids} > sample_ids.txt
        echo ~{sep=' ' family_ids} > family_ids.txt
        echo ~{sep=' ' gvcf_paths} > gvcf_paths.txt
        echo ~{sep=' ' gvcf_index_paths} > gvcf_index_paths.txt
        echo ~{sep=' ' pedigrees} > pedigrees.txt
        echo ~{sep=' ' reported_sexes} > reported_sexes.txt
        echo ~{family_id} > family_id.txt

        python3 <<CODE
        import os
        import shutil
        import subprocess


        def process_family(sample_ids, family_ids, gvcf_paths, family_id):
            print(f"Processing family: {family_id}.")
            
            family_gvcfs = []

            for sample_id, fam_id, gvcf_path in zip(sample_ids, family_ids, gvcf_paths):
                if fam_id == family_id:
                    family_gvcfs.append(gvcf_path)
                    print(f"Added gVCF for sample {sample_id}: {gvcf_path}")

                    # Index the gVCF file using bcftools
                    try:
                        subprocess.run(['bcftools', 'index', '-t', gvcf_path], check=True)
                        print(f"Indexed gVCF for sample {sample_id}: {gvcf_path}")
                    except subprocess.CalledProcessError as e:
                        print(f"Error indexing gVCF for sample {sample_id}: {e}")

            print(f"There are {len(family_gvcfs)} samples in this family.")

            return family_gvcfs

        def write_family_info(sample_ids, pedigrees, reported_sexes, family_ids, family_id):
            filename = f"family_info_{family_id}.txt"
            header = "sample_id\tpedigree\treported_sex\tsidr_family_id\n"

            with open(filename, "w") as f:
                # Write the header
                f.write(header)              
                # Iterate over the lists and filter by family_id
                for sample_id, pedigree, reported_sex, fam_id in zip(sample_ids, pedigrees, reported_sexes, family_ids):
                    if fam_id == family_id:
                        line = f"{sample_id}\t{pedigree}\t{reported_sex}\t{family_id}\n"
                        f.write(line)

        with open("sample_ids.txt", "r") as f:
                sample_ids = f.read().strip().split(' ')
        with open("family_ids.txt", "r") as f:
                family_ids = f.read().strip().split(' ')
        with open("gvcf_paths.txt", "r") as f:
                gvcf_paths = f.read().strip().split(' ')
        with open("pedigrees.txt", "r") as f:
                pedigrees = f.read().strip().split(' ')
        with open("reported_sexes.txt", "r") as f:
                reported_sexes = f.read().strip().split(' ')    
        with open("family_id.txt", "r") as f:
                family_id = f.read().strip()
        

        family_gvcfs = process_family(sample_ids, family_ids, gvcf_paths, family_id)
        write_family_info(sample_ids, pedigrees, reported_sexes, family_ids, family_id)

        family_gvcfs_str = ' '.join(family_gvcfs)

        vcf_merge_command = f"vcf-merge {family_gvcfs_str} > merged_family_{family_id}.vcf"
        bgzip_command = f"bgzip merged_family_{family_id}.vcf"
        bcftools_index_command = f"bcftools index -t merged_family_{family_id}.vcf.gz"

        subprocess.run(vcf_merge_command, shell=True, check=True)
        subprocess.run(bgzip_command, shell=True, check=True)
        subprocess.run(bcftools_index_command, shell=True, check=True)

        CODE
        >>>

    output {
        File merged_vcf = "merged_family_~{family_id}.vcf.gz"
        File merged_vcf_index = "merged_family_~{family_id}.vcf.gz.tbi"
        File family_info = "family_info_~{family_id}.txt"
    }

    runtime {
        docker: "us.gcr.io/tag-public/peddy-analysis:v1"
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}



task UpdateFamFile {
    input {
        File fam_file
        String family_id
        File known_trio_info
        Int memory = 16
        Int disk_size = 16
    }

    command <<<
        python3 <<CODE

        import pandas as pd

        def update_fam_file(fam_file_path, known_trio_info_path, output_fam_file_path):
            # Read the .fam file into a DataFrame
            fam_file_df = pd.read_csv(fam_file_path, sep=' ', header=None, names=['FamilyID', 'SampleID', 'FatherID', 'MotherID', 'Sex', 'Phenotype'])

            known_trio_info_df = pd.read_csv(known_trio_info_path, sep='\t')
            sample_info = known_trio_info_df.set_index('sample_id').to_dict(orient='index')

            # Ensure the correct data type for columns
            fam_file_df[['FamilyID', 'FatherID', 'MotherID']] = fam_file_df[['FamilyID', 'FatherID', 'MotherID']].astype(str)

            # Update the .fam file DataFrame based on the trio information
            for i, row in fam_file_df.iterrows():
                sample_id = row['SampleID']
                if sample_id in sample_info:
                    info = sample_info[sample_id]
                    fam_file_df.at[i, 'FamilyID'] = info['sidr_family_id']
                    if info['pedigree'] == 'Proband':
                        father_id = known_trio_info_df[(known_trio_info_df['sidr_family_id'] == info['sidr_family_id']) & (known_trio_info_df['pedigree'] == 'Father')]['sample_id'].values[0]
                        mother_id = known_trio_info_df[(known_trio_info_df['sidr_family_id'] == info['sidr_family_id']) & (known_trio_info_df['pedigree'] == 'Mother')]['sample_id'].values[0]
                        fam_file_df.at[i, 'FatherID'] = father_id
                        fam_file_df.at[i, 'MotherID'] = mother_id
                    else:
                        fam_file_df.at[i, 'FatherID'] = '0'
                        fam_file_df.at[i, 'MotherID'] = '0'
                    fam_file_df.at[i, 'Sex'] = 1 if info['reported_sex'] == 'Male' else 2

            # Save the updated .fam file
            fam_file_df.to_csv(output_fam_file_path, sep=' ', header=False, index=False)

        # Update and Save the updated .fam file
        update_fam_file("~{fam_file}", "~{known_trio_info}", "updated_~{family_id}.fam")
        CODE
    >>>

    output {
        File updated_fam_file = "updated_~{family_id}.fam"
    }

    runtime {
        docker: "us.gcr.io/tag-public/peddy-analysis:v1"
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task RunPlink {
    input {
        String family_id
        File merged_vcf
        File merged_vcf_index
        Int memory = 16
        Int disk_size = 16
    }

    command <<<
        # Run PLINK on the merged VCF file
        plink --vcf ~{merged_vcf} --make-bed --out ~{family_id} --allow-extra-chr
    >>>

    output {
        File binary_bim = "~{family_id}.bim"
        File binary_fam = "~{family_id}.fam"
        File binary_bed = "~{family_id}.bed"
    }

    runtime {
        docker: "us.gcr.io/tag-public/peddy-analysis:v1"
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task RunPeddy {
    input {
        String prefix
        File merged_vcf
        File merged_vcf_index
        File fam_file
        String reference_genome = "hg38"
        Int memory = 16
        Int disk_size = 16
    }
    command <<<        
        peddy -p 4 --plot --prefix ~{prefix} ~{merged_vcf} ~{fam_file} --sites ~{reference_genome}

    >>>
    output {
        File ancestry_assignment_plot = "~{prefix}.pca_check.png"
        File pedigree_prediction_stats = "~{prefix}.ped_check.csv"
        File pedigree_prediction_plot = "~{prefix}.ped_check.png"
    }

    runtime {
        docker: "us.gcr.io/tag-public/peddy-analysis:v1"
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task MergePeddyResults {
    input {
        Array[File] peddy_results
        Int memory = 16
        Int disk_size = 16
    }
    command <<<
        # Concatenate all the Peddy results
        echo "Merging all Peddy outputs into one"
        head -n 1 ~{peddy_results[0]} > merged_peddy_results.csv

        for file in ~{sep=' ' peddy_results}; do
            tail -n +2 $file >> merged_peddy_results.csv
        done
    
    >>>

    output {
        File merged_peddy_results = "merged_peddy_results.csv"
    }

    runtime {
        docker: "us.gcr.io/tag-public/peddy-analysis:v1"
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task PlotPeddyResults {
    input {
        File merged_peddy_results
        Int memory = 8
        Int disk_size = 16
    }

    command <<<
        # Create a Python script to generate the plot
        python3 <<CODE
        import pandas as pd
        import matplotlib.pyplot as plt

        def plot_peddy_results(input_file, output_file):
            df = pd.read_csv(input_file)
            fig, ax = plt.subplots(figsize=(15, 12))

            # Define color mapping
            colors = {'True': 'red', 'False': 'green'}

            for idx, row in df.iterrows():
                color = colors.get(str(row['parent_error']), 'black')
                ax.scatter(row['pedigree_relatedness'], row['rel'], color=color, marker='o')
                if str(row['parent_error']) == 'True':
                    ax.text(row['pedigree_relatedness'], row['rel'], f"{row['sample_a']}-{row['sample_b']}", fontsize=9)

            plt.title('Peddy Pedigree Predictions')
            plt.xlabel('Pedigree Relatedness')
            plt.ylabel('Peddy Infered Relatedness')
            plt.xlim(0, 1)
            plt.ylim(0, 1)
            plt.plot([0, 1], [0, 1], color='gray', linestyle='--')  # Add y = x line
            plt.savefig(output_file, dpi=300)

        plot_peddy_results("~{merged_peddy_results}", "peddy_prediction_plot.png")
        CODE

    >>>

    output {
        File  all_family_peddy_prediction_plot = "peddy_prediction_plot.png"
    }

    runtime {
        docker: "us.gcr.io/tag-public/peddy-analysis:v1"
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}