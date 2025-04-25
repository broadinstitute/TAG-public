version 1.0


workflow Peddy_AnalyzeFamilialRelatedness {
    input {
        Array[String] sample_ids
        Array[String] family_ids
        Array[String] gvcf_paths
        Array[String] gvcf_index_paths
        Array[String] pedigrees
        Array[String] reported_sexes
        File interval_list
        File reference_fasta
        File reference_fasta_index
        File reference_dict
    }
        

    call AnalyzeAndCheckFamilySamples {
        input:  
        sample_ids = sample_ids, 
        family_ids = family_ids
    }
    call FilterSingleSampleFamilies {
        input:
            unique_family_ids = AnalyzeAndCheckFamilySamples.unique_family_ids,
            sample_ids = sample_ids,
            family_ids = family_ids,
            gvcfs = gvcf_paths,
            gvcf_indexes = gvcf_index_paths,
            pedigrees = pedigrees,
            reported_sexes = reported_sexes
    }
    scatter (index in range(length(FilterSingleSampleFamilies.passing_family_ids))){
        call GroupFamilyGVCFs {
            input: 
                grouped_per_family_gvcf = FilterSingleSampleFamilies.grouped_per_family_gvcf[index]
        }
        call ProcessFamilyGVCFs {
            input:
                gvcf_paths = GroupFamilyGVCFs.gvcf_paths,
                family_id = GroupFamilyGVCFs.family_id
        }
        call CombineGVCFs {
            input:
                output_prefix = GroupFamilyGVCFs.family_id,
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                reference_dict = reference_dict,
                family_gvcfs = ProcessFamilyGVCFs.family_gvcfs,
                interval_list = interval_list

        }
        call GenotypeGVCFs {
            input: 
                combined_gvcf = CombineGVCFs.combined_gvcf,
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                reference_dict = reference_dict,
                output_prefix = GroupFamilyGVCFs.family_id

        }
        call RunPlink {
            input:
                family_id = GroupFamilyGVCFs.family_id,
                merged_gvcf = GenotypeGVCFs.genotyped_gvcf,
                merged_gvcf_index = GenotypeGVCFs.genotyped_gvcf_index
        }
        call UpdateFamFile {
            input:
                fam_file = RunPlink.binary_fam,
                known_trio_info = FilterSingleSampleFamilies.family_info[index],
                family_id = GroupFamilyGVCFs.family_id
        }
        call RunPeddy {
            input:
                prefix = GroupFamilyGVCFs.family_id,
                merged_gvcf = GenotypeGVCFs.genotyped_gvcf,
                merged_gvcf_index = GenotypeGVCFs.genotyped_gvcf_index,
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



    output {
        File family_composition_log =  AnalyzeAndCheckFamilySamples.family_composition_log
        Int num_single_sample_families =  AnalyzeAndCheckFamilySamples.num_single_sample_families
        Array[String] processed_families = FilterSingleSampleFamilies.passing_family_ids
        Array[String] unprocessed_family_ids = FilterSingleSampleFamilies.unprocessed_family_ids
        Array[File] merged_gvcf_files = GenotypeGVCFs.genotyped_gvcf
        Array[File] merged_gvcf_indices = GenotypeGVCFs.genotyped_gvcf_index
        Array[File] family_info_files = FilterSingleSampleFamilies.family_info
        Array[File] updated_fam_file = UpdateFamFile.updated_fam_file
        Float concordance = MergePeddyResults.concordance
        File merged_peddy_results = MergePeddyResults.merged_peddy_results
        File all_family_peddy_prediction_plot = PlotPeddyResults.all_family_peddy_prediction_plot
    }
    
}

task AnalyzeAndCheckFamilySamples {
    input {
        Array[String] sample_ids
        Array[String] family_ids
        Int memory = 16
        Int disk_size = 16
        String? docker_override
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
            num_large_families = (family_counts >= 4).sum() 

            with open(log_file, 'w') as log:
                # Check for families with only one sample
                single_sample_families = family_counts[family_counts == 1]
                if not single_sample_families.empty:
                    for family_id in single_sample_families.index:
                        sample_id = df[df['sidr_family_id'] == family_id]['entity:sample_id'].values[0]
                        log.write(f'Warning: Family {family_id} has only one sample: {sample_id}\n')

                log.write(f'This sample set has {num_trios} Trios and {num_duos} Duos.\n')
                log.write(f'The number of families with size >= 4 is {num_large_families}.\n') 
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
        docker: select_first([docker_override, "us.gcr.io/tag-public/peddy-analysis:v1"])
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task FilterSingleSampleFamilies {
    input {
        Array[String] unique_family_ids
        Array[String] sample_ids
        Array[String] family_ids
        Array[String] gvcfs
        Array[String] gvcf_indexes
        Array[String] pedigrees
        Array[String] reported_sexes
        String? docker_override
    }

    command <<<
        echo ~{sep=' ' sample_ids} > sample_ids.txt
        echo ~{sep=' ' unique_family_ids} > unique_family_ids.txt
        echo ~{sep=' ' family_ids} > family_ids.txt
        echo ~{sep=' ' gvcfs} > gvcf_paths.txt
        echo ~{sep=' ' gvcf_indexes} > gvcf_index_paths.txt
        echo ~{sep=',' pedigrees} > pedigrees.txt
        echo ~{sep=' ' reported_sexes} > reported_sexes.txt

        # Filter out single-sample families and generate the required files
        python3 <<CODE

        import os

        with open("sample_ids.txt", "r") as f:
            sample_ids = f.read().strip().split(' ')
        with open("unique_family_ids.txt", "r") as f:
            unique_family_ids = f.read().strip().split(' ')
        with open("family_ids.txt", "r") as f:
            family_ids = f.read().strip().split(' ')
        with open("gvcf_paths.txt", "r") as f:
            gvcfs = f.read().strip().split(' ')
        with open("gvcf_index_paths.txt", "r") as f:
            gvcf_indexes = f.read().strip().split(' ')
        with open("pedigrees.txt", "r") as f:
            pedigrees = f.read().strip().split(',')
        with open("reported_sexes.txt", "r") as f:
            reported_sexes = f.read().strip().split(' ')

        pass_family_ids = []
        excluded_family_ids = []

        for family_id in unique_family_ids:
            family_sample_ids = [sample_id for sample_id, fid in zip(sample_ids, family_ids) if fid == family_id]
            family_gvcfs = [gvcf for gvcf, fid in zip(gvcfs, family_ids) if fid == family_id]
            family_gvcf_indexes = [gvcf_index for gvcf_index, fid in zip(gvcf_indexes, family_ids) if fid == family_id]
            family_pedigrees = [pedigree for pedigree, fid in zip(pedigrees, family_ids) if fid == family_id]
            family_reported_sexes = [reported_sex for reported_sex, fid in zip(reported_sexes, family_ids) if fid == family_id]

            if len(family_sample_ids) > 1:
                pass_family_ids.append(family_id)

                with open(f"family_info_{family_id}.txt", "w") as family_info_file, \
                     open(f"grouped_per_family_gvcf_{family_id}.txt", "w") as grouped_gvcf_file:
                    family_info_file.write("sample_id\tpedigree\treported_sex\tsidr_family_id\n")
                    grouped_gvcf_file.write("family_id\tsample_id\tgvcf_path\tgvcf_index_path\n")

                    for sample_id, gvcf, gvcf_index, pedigree, reported_sex in zip(family_sample_ids, family_gvcfs, family_gvcf_indexes, family_pedigrees, family_reported_sexes):
                        # Write to individual family_info.txt
                        family_info_file.write(f"{sample_id}\t{pedigree}\t{reported_sex}\t{family_id}\n")
                        # Write to individual grouped_per_family_gvcf.txt
                        grouped_gvcf_file.write(f"{family_id}\t{sample_id}\t{gvcf}\t{gvcf_index}\n")
            else:
                excluded_family_ids.append(family_id)

        # Write the excluded and passing family IDs
        with open("excluded_family_ids.txt", "w") as f:
            f.write("\n".join(excluded_family_ids))

        with open("passing_family_ids.txt", "w") as f:
            f.write("\n".join(pass_family_ids))

        CODE
    >>>

    output {
        Array[File] family_info = glob("family_info_*.txt")
        Array[File] grouped_per_family_gvcf = glob("grouped_per_family_gvcf_*.txt")
        Array[String] unprocessed_family_ids = read_lines("excluded_family_ids.txt")
        Array[String] passing_family_ids = read_lines("passing_family_ids.txt")
    }

    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/peddy-analysis:v1"])
        memory: "8 GB"
        disks: "local-disk 16 HDD"
    }
}

task GroupFamilyGVCFs {
    input {
        File grouped_per_family_gvcf
        Int disk_size = 10
        Int memory = 32
        Int? buffer_disk_size
        String? docker_override
    }

    command <<<
        python3 <<CODE

        import os

        family_file = "~{grouped_per_family_gvcf}"

        sample_ids = []
        gvcf_paths = []
        gvcf_index_paths = []
        consistent_family_id = None
        family_id_set = set()

        # Parse the grouped gVCF file
        with open(family_file, "r") as f:
            next(f) 
            for line in f:
                fields = line.strip().split("\t")
                family_id = fields[0]
                family_id_set.add(family_id)
                sample_ids.append(fields[1])
                gvcf_paths.append(fields[2])
                gvcf_index_paths.append(fields[3])

        # Ensure consistent family_id in the file
        if len(family_id_set) == 1:
            consistent_family_id = family_id_set.pop()
        else:
            raise ValueError("Inconsistent family_id values found in the file.")

        # Save parsed data to intermediate files
        with open("consistent_family_id.txt", "w") as f:
            f.write(consistent_family_id)
        with open("gvcf_paths.txt", "w") as f:
            f.write("\n".join(gvcf_paths))
        with open("gvcf_index_paths.txt", "w") as f:
            f.write("\n".join(gvcf_index_paths))
        with open("sample_ids.txt", "w") as f:
            f.write("\n".join(sample_ids))

        CODE
    >>>

    output {
        String family_id = read_string("consistent_family_id.txt")
        Array[File] gvcf_paths = read_lines("gvcf_paths.txt")
        Array[File] gvcf_index_paths = read_lines("gvcf_index_paths.txt")
        Array[String] sample_ids = read_lines("sample_ids.txt")
    }

    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/peddy-analysis:v1"])
        memory: memory + " GB"
        disks: "local-disk " + (disk_size + select_first([buffer_disk_size, 0])) + " SSD"
    }
}

task ProcessFamilyGVCFs {
    input {
        Array[File] gvcf_paths
        String family_id
        Int disk_size = 32
        Int memory = 32
        Int? buffer_disk_size
        String? docker_override
    }

    command <<<
        python3 <<CODE

        import os
        import subprocess
        from urllib.parse import urlparse

        # read in gvcf_paths as input
        gvcf_paths = "~{sep=' ' gvcf_paths}".split(" ")

        # Process and index gVCF files for the family
        def process_family(gvcf_paths, family_id):
            print(f"Processing family: {family_id}.")
            
            family_gvcfs = []
            
            for gvcf_path in gvcf_paths:
                gvcf_base_name = os.path.basename(urlparse(gvcf_path).path)
                filtered_gVCF_path = f"filtered_{gvcf_base_name}"
                print(f"Starting to filter: {gvcf_path}")
                try:
                    # Filter the gVCF file and ensure the ALT column is all <NON_REF>
                    import subprocess

                    # Filter the gVCF file and ensure the ALT column is all <NON_REF>
                    subprocess.run(f'bcftools view -i \'ALT=="<NON_REF>"\' {gvcf_path} -Oz -o {filtered_gVCF_path}', shell=True, check=True)
                    family_gvcfs.append(filtered_gVCF_path)

                except subprocess.CalledProcessError as e:
                    print(f"Error filtering gVCF for {gvcf_path}: {e.stderr}")

            print(f"Processed {len(family_gvcfs)} samples in this family.")
            return family_gvcfs


        # Process per family gVCFs
        family_gvcfs = process_family(gvcf_paths, "~{family_id}")


        CODE
    >>>

    output {
        Array[File] family_gvcfs = glob("filtered_*.gvcf.gz")
    }

    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/peddy-analysis:v1"])
        memory: memory + " GB"
        disks: "local-disk " + (disk_size + select_first([buffer_disk_size, 0])) + " SSD"
    }
}


task CombineGVCFs {
    input {
        String output_prefix
        File reference_fasta
        File reference_fasta_index
        File reference_dict 
        File interval_list
        Array[File] family_gvcfs
        Int memory = 64
        Int disk_size = ceil(size(reference_fasta, 'GB') +
                                size(reference_fasta_index, 'GB') +
                                size(reference_dict, 'GB') +
                                (length(family_gvcfs) * 5)) + disk_pad
        Int disk_pad = 0
        Int preemptible = 2
    }


    command <<<

        # Index gVCF files using bcftools
        for gvcf in ~{sep=' ' family_gvcfs}; do
            echo "Indexing $gvcf"
            bcftools index -t $gvcf
        done       

        gatk CombineGVCFs \
        -R ~{reference_fasta} \
        ~{sep=' ' prefix("--variant ", family_gvcfs)} \
        -O ~{output_prefix}.g.vcf.gz \
        -L ~{interval_list}

    >>>

    output {
      File combined_gvcf = "~{output_prefix}.g.vcf.gz"
    }
    
    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.6.0.0"
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
        maxRetries: 2
        preemptible: "${preemptible}"
    }
}

task GenotypeGVCFs {
    input {
        File combined_gvcf
        File reference_fasta
        File reference_fasta_index
        File reference_dict 
        String output_prefix
        Int disk_pad = 0
        Int disk_size = ceil(size(combined_gvcf, "GB")) + disk_pad
        Int memory = 32
        Int preemptible = 2
    }


    command <<<
    
     bcftools index -t ~{combined_gvcf} > ~{output_prefix}.g.vcf.gz.tbi
     
     # genotype grouped a multiple-sample gVCF
     gatk GenotypeGVCFs \
       -R ~{reference_fasta} \
       -V ~{combined_gvcf} \
       --read-index ~{output_prefix}.g.vcf.gz.tbi \
       -O ~{output_prefix}.genotyped.vcf.gz
            
    >>>

    output {
        File genotyped_gvcf = "~{output_prefix}.genotyped.vcf.gz"
        File genotyped_gvcf_index = "~{output_prefix}.genotyped.vcf.gz.tbi"
    }

    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.6.0.0"
        memory: memory + "GB"
        cpu: "2"
        disks: "local-disk " + disk_size + " HDD"
        maxRetries: 2
        preemptible: "${preemptible}"
    }
}



task RunPlink {
    input {
        String family_id
        File merged_gvcf
        File merged_gvcf_index
        Int memory = 16
        Int disk_size = 16
        String? docker_override
        
    }

    command <<<
        # Run PLINK on the merged VCF file
        plink --vcf ~{merged_gvcf} --make-bed --out ~{family_id} --allow-extra-chr
    >>>

    output {
        File binary_bim = "~{family_id}.bim"
        File binary_fam = "~{family_id}.fam"
        File binary_bed = "~{family_id}.bed"
    }

    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/peddy-analysis:v1"])
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
        String? docker_override

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
                        # father_id_row = known_trio_info_df[(known_trio_info_df['sidr_family_id'] == info['sidr_family_id']) & (known_trio_info_df['pedigree'] == 'Father')]
                        # mother_id_row = known_trio_info_df[(known_trio_info_df['sidr_family_id'] == info['sidr_family_id']) & (known_trio_info_df['pedigree'] == 'Mother')]
                        # This is to handle cases like "Biological Father/Mother"
                        father_id_row = known_trio_info_df[
                            (known_trio_info_df['sidr_family_id'] == info['sidr_family_id']) &
                            (known_trio_info_df['pedigree'].str.contains('Father', case=False, na=False))
                        ]
                        mother_id_row = known_trio_info_df[
                            (known_trio_info_df['sidr_family_id'] == info['sidr_family_id']) &
                            (known_trio_info_df['pedigree'].str.contains('Mother', case=False, na=False))
                        ]
                        # This is to handle if only one parent present
                        father_id = '0'
                        mother_id = '0'                
                        if not father_id_row.empty:
                            father_id = father_id_row['sample_id'].values[0]
                        if not mother_id_row.empty:
                            mother_id = mother_id_row['sample_id'].values[0]
                        
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
        docker: select_first([docker_override, "us.gcr.io/tag-public/peddy-analysis:v1"])
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}



task RunPeddy {
    input {
        String prefix
        File merged_gvcf
        File merged_gvcf_index
        File fam_file
        String? reference_genome
        Int memory = 16
        Int disk_size = 16
        String? docker_override
    }
    command <<<        
        peddy -p 4 --plot --prefix ~{prefix} ~{merged_gvcf} ~{fam_file}  ~{if defined(reference_genome) then "--sites " + reference_genome else ""}

    >>>
    output {
        File ancestry_assignment_plot = "~{prefix}.pca_check.png"
        File pedigree_prediction_stats = "~{prefix}.ped_check.csv"
        File pedigree_prediction_plot = "~{prefix}.ped_check.png"
    }

    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/peddy-analysis:v1"])
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task MergePeddyResults {
    input {
        Array[File] peddy_results
        Int memory = 16
        Int disk_size = 16
        String? docker_override
    }
    command <<<
        # Concatenate all the Peddy results
        echo "Merging all Peddy outputs into one"
        head -n 1 ~{peddy_results[0]} > merged_peddy_results.csv

        for file in ~{sep=' ' peddy_results}; do
            tail -n +2 $file >> merged_peddy_results.csv
        done


        # Calculate the prediction concordance
        python3 <<CODE
        import pandas as pd

        df = pd.read_csv("merged_peddy_results.csv")

        # Calculate the concordance as the ratio of "False" in the parent_error column
        total_rows = len(df)
        false_count = df['parent_error'].astype(str).value_counts().get('False', 0)
        concordance = false_count / total_rows if total_rows > 0 else 0
        concordance = round(concordance, 2) # round by 2 decimals

        with open("concordance.txt", "w") as f:
            f.write(str(concordance))
        CODE
    
    >>>

    output {
        File merged_peddy_results = "merged_peddy_results.csv"
        Float concordance = read_float("concordance.txt")
    }

    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/peddy-analysis:v1"])
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task PlotPeddyResults {
    input {
        File merged_peddy_results
        Int memory = 8
        Int disk_size = 16
        String? docker_override
    }

    command <<<

        python3 <<CODE

        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt


        df = pd.read_csv("~{merged_peddy_results}")

        df['concordance'] = df['parent_error'].apply(lambda x: 'Disconcordant' if x else 'Concordant')

        custom_palette = {'Concordant': 'green', 'Disconcordant': 'red'}
        total_num_points = len(df)
        mean_n = df['n'].mean()
        unique_concordance = df['concordance'].unique()
        g = sns.catplot(x="concordance", y="rel_difference", kind="swarm", data=df, height=6, aspect=1.5,  palette=custom_palette)

        for i, label in enumerate(unique_concordance):
            num_points = len(df[df['concordance'] == label])
            percentage = (num_points / len(df)) * 100
            plt.text(i, -0.9, f'{num_points} ({percentage:.2f}%)', 
                    ha='center', fontsize=12, color='gray')


        plt.ylim(-1, 1)

        plt.title('Peddy Concordance with Relatedness Difference', y=1)
        plt.xlabel('Prediction VS Known')
        plt.ylabel('Rel Difference')
        plt.tight_layout()
        plt.axhline(y=0, color='grey', linestyle='--')

        # Add the mean variants count as a label
        plt.suptitle(f'( # of sample pairs = {total_num_points} | Mean Variants (n) = {mean_n:.2f} )', fontsize=12, ha='center', color='black', y=0.94)

        plt.savefig("peddy_prediction_plot.png", dpi = 300)

        CODE

    >>>

    output {
        File  all_family_peddy_prediction_plot = "peddy_prediction_plot.png"
    }

    runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/peddy-analysis:v1"])
        memory: memory + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}