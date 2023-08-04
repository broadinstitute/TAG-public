version 1.0
workflow AnnotateBed{
    input {
        File annotate_script
        File gencode_annotation
        File bed_to_annotate
        String output_prefix
        File gene_bed
        Boolean generate_gene_base_count
        File? gene_list
        File gene_base_count_script
        Int diskGB = 50
        Array[File] bam_files
        Array[File] bam_indices
        File bed_file
        Array[File] sample_ids
        Int memory_gb
        Int disk_size
    }
    call GenerateAnnotation {
        input:
        script = annotate_script,
        gencode_annotation = gencode_annotation, 
        bed_to_annotate = bed_to_annotate,
        output_prefix = output_prefix,
        gene_list = gene_list,
        gene_bed = gene_bed
    }
    if ((generate_gene_base_count) && defined(gene_base_count_script) && defined(gene_bed)) {
        call CountGeneBases {
            input:
            script = gene_base_count_script,
            bed_to_annotate = bed_to_annotate,
            gene_bed = gene_bed,
            gene_list = gene_list,
            grouped_by_gene = GenerateAnnotation.grouped_by_gene,
            output_prefix = output_prefix
            }
        }
    call samtools_coverage {
        input:
        bam_files = bam_files,
        bam_indices = bam_indices,
        sample_ids = sample_ids,
        bed_file = bed_to_annotate,
        memory_gb = memory_gb,
        disk_size = disk_size
    }
    call summarize_coverage {
        input:
        coverage_output = samtools_coverage.coverage_output
    }
    call generate_clinvar_results {
        input:
        bed_to_annotate = bed_to_annotate,
        memory_gb = memory_gb,
        disk_size = disk_size

    }
    call concatenate_results {
        input:
        clinvar_annotation = generate_clinvar_results.clinvar_annotation, 
        samtools_coverage_file = summarize_coverage.samtools_coverage_summary,
        annotation_file = GenerateAnnotation.annotation_per_interval,
        output_prefix = output_prefix, 
        memory_gb = memory_gb,
        disk_size = disk_size
    }
    output {
        File ungrouped_annotation = GenerateAnnotation.ungrouped_annotation
        File annotation_per_interval = GenerateAnnotation.annotation_per_interval
        Int intergenic_base_count = read_int(GenerateAnnotation.intergenic_base_count_file)
        Int coding_base_count = read_int(GenerateAnnotation.coding_base_count_file)
        File grouped_by_gene = GenerateAnnotation.grouped_by_gene 
        Array[Int?] genes_involved = CountGeneBases.genes_involved
        Array[File?] gene_base_count = CountGeneBases.gene_base_count
        File samtools_coverage_summary = summarize_coverage.samtools_coverage_summary
        File clinvar_annotation = generate_clinvar_results.clinvar_annotation
        File integrated_annotation_file = concatenate_results.integrated_annotation_file
    }
    
}
task GenerateAnnotation {
    input {
        File script
        File gencode_annotation
        File bed_to_annotate
        File gene_bed
        File? gene_list
        String output_prefix
        Int maxRetries = 1
        Int preemptible = 3
        Int diskGB = 50
    }
    command {
        python3 ~{script} --annotation ~{gencode_annotation} --bed ~{bed_to_annotate} --gene_bed ~{gene_bed} ~{'--gene_list ' + gene_list} --output_prefix ~{output_prefix}
    }
    runtime {
        docker: "us.gcr.io/tag-team-160914/annotatebed:test"
        preemptible: preemptible
        disks: "local-disk ~{diskGB} HDD"
        memory: "4GB"
        maxRetries: maxRetries
    }
    output {
        File ungrouped_annotation = "~{output_prefix}.ungrouped.annotated.txt"
        File annotation_per_interval = "~{output_prefix}.grouped_by_interval.annotated.txt"
        File grouped_by_gene = "~{output_prefix}.grouped_by_gene.txt"
        File coding_base_count_file= "coding_base_count.txt"
        File intergenic_base_count_file = "intergenic_base_count.txt"
    }
}

task CountGeneBases {
    input{
        File script
        File gene_bed
        File bed_to_annotate
        File? gene_list
        File grouped_by_gene = grouped_by_gene
        Int maxRetries = 1
        Int preemptible = 3
        Int diskGB = 50
        String output_prefix
    }
    command {
            python3 ~{script} --gene_bed ~{gene_bed} --bed ~{bed_to_annotate} --grouped_by_gene ~{grouped_by_gene} ~{'--gene_list ' + gene_list} --output_prefix ~{output_prefix}
    }
    runtime {
        docker: "us.gcr.io/tag-team-160914/annotatebed:test"
        preemptible: preemptible
        disks: "local-disk ~{diskGB} HDD"
        memory: "4GB"
        maxRetries: maxRetries
    }
    output {
        Int genes_involved = read_int(stdout())
        File? gene_base_count = "~{output_prefix}.gene_base_count.txt"
    }
}
task samtools_coverage {
    input {
        Array[File] bam_files
        Array[File] bam_indices
        File bed_file
        Array[File] sample_ids
        Int memory_gb
        Int disk_size
    }
    command <<<
    for ((i = 0; i < length(bam_files); i++))
    do
        sample_id = sample_ids[i]
        bam_file = bam_files[i]
        tmp_file = "tmp_${sample_id}"

        for region in $(awk '{print $1":"$2"-"$3}' ~{bed_file})
        do
            samtools coverage $bam_file -r ${region} >> $tmp_file
        done

        awk '/startpos/ && c++ > 0 {next} 1' $tmp_file > "${sample_id}.coverage.txt"
    done
>>>
    output {
        Array[File] coverage_output = glob("*.coverage.txt")
    }
    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.5.7-2021-06-09_16-47-48Z"
        memory: memory_gb + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
task summarize_coverage {
    input{
        Array[File] coverage_output
        Int memory_gb
        Int disk_size
    }
    command {
    python3 <<CODE
import pandas as pd
file_paths = ${sep=' ' coverage_output}
counts_dict = {}
columns_to_count = ['meandepth', 'meanbaseq', 'meanmapq']
dfs = [pd.read_csv(file_path, sep='\t', usecols=[
                   '#rname', 'startpos', 'endpos'] + columns_to_count) for file_path in file_paths]
# looping over each line and each file
for df in dfs:
    for index, row in df.iterrows():
        rname = row['#rname']
        startpos = row['startpos']
        endpos = row['endpos']
        # Check if the condition (column > 20) is satisfied for this row in this DataFrame for each column
        for column in columns_to_count:
            column_gt_20 = row[column] > 20
            if column_gt_20:
                key = (rname, startpos, endpos)
                counts_dict.setdefault(key, {}).setdefault(column, 0)
                counts_dict[key][column] += 1

results_df = pd.DataFrame.from_dict(counts_dict, orient='index').reset_index()
results_df.rename(columns={'level_0': 'chrom', 'level_1': 'start', 'level_2': 'end', 'meandepth': 'count_meandepth_gt_20',
                  'meanbaseq': 'count_meanbaseq_gt_20', 'meanmapq': 'count_meanmapq_gt_20'}, inplace=True)
print(results_df)

results_df['chrom'] = results_df['chrom'].str.replace(
    r'^chr', '', regex=True).astype(int)

# Sort the DataFrame by genomic coordinates (first column, then second column)
results_df = results_df.sort_values(by=['chrom', 'start'])

results_df = results_df.reset_index(drop=True)
results_df['chrom'] = results_df['chrom'].apply(lambda x: f'chr{x}')
results_df['count_meandepth_gt_20'] = results_df['count_meandepth_gt_20'].fillna(
    0).astype(int)
results_df['count_meanbaseq_gt_20'] = results_df['count_meanbaseq_gt_20'].fillna(
    0).astype(int)
results_df['count_meanmapq_gt_20'] = results_df['count_meanmapq_gt_20'].fillna(
    0).astype(int)

# This code below is to figure out the failing reasons with criteria: 80-20-20-20
num_files = len(file_paths)
results_df['%_of_samples_with_low_depth'] = ((1 -
                                             (results_df['count_meandepth_gt_20'] / num_files)) * 100).astype(int)
results_df['%_of_samples_with_low_baseq'] = ((1 -
                                             (results_df['count_meanbaseq_gt_20'] / num_files)) * 100).astype(int)
results_df['%_of_samples_with_low_mapq'] = ((1 -
                                            (results_df['count_meanmapq_gt_20'] / num_files)) * 100).astype(int)

results_df['failing_reason'] = ''
results_df.loc[results_df['%_of_samples_with_low_depth']
               > 80, 'failing_reason'] += 'depth,'
results_df.loc[results_df['%_of_samples_with_low_baseq']
               > 80, 'failing_reason'] += 'baseq,'
results_df.loc[results_df['%_of_samples_with_low_mapq']
               > 80, 'failing_reason'] += 'mapq,'
results_df['failing_reason'] = results_df['failing_reason'].str.rstrip(',')
results_df.to_csv('samtools_coverage_summary.txt', index=False)

    }
    output{
    File samtools_coverage_summary = 'samtools_coverage_summary.txt'
    }
    runtime {
    	docker: "us.gcr.io/tag-team-160914/integrate_annotation:test"
        memory: memory_gb + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
task generate_clinvar_results{
    input {
        File bed_to_annotate
        Int memory_gb
        Int disk_size
    }
    command <<<

        bcftools query -f '%CHROM\t%POS\t%INFO/ALLELEID\t%INFO/CLNHGVS\t%INFO/CLNREVSTAT\t%INFO/CLNSIG\n' -i 'INFO/CLNSIG="Likely_pathogenic" || INFO/CLNSIG="Pathogenic" || INFO/CLNSIG="Pathogenic/Likely_pathogenic"' -R $bed_to_annotate clinvar.vcf.gz >  clinvar_annotation.txt
    >>>
    output {
        File clinvar_annotation = "clinvar_annotation.txt"
    }
    runtime {
        docker: "us.gcr.io/tag-team-160914/integrate_annotation:test"
        memory: memory_gb + "GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
task concatenate_results {
    input{
        File clinvar_annotation
        File samtools_coverage_file
        File annotation_file
        String output_prefix
        Int memory_gb
        Int disk_size
    }
    command {
        python3 /scripts/aggregation_script.py --clinvar_file ~{clinvar_annotation} --samtools_coverage_file ~{samtools_coverage_file} --annotation_file ~{annotation_file} --output_prefix ~{output_prefix}
    }
    output {
        File integrated_annotation_file = "~{output_prefix}.integrated_annotation.txt"
    }
    runtime {
        docker: "us.gcr.io/tag-team-160914/integrate_annotation:test"
        memory: memory_gb + "GB"
        disks: "local-disk " + disk_size + " HDD"

    }
}