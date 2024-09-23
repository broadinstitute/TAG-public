import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description='Add input files.')
    parser.add_argument('--clinvar_file', type=str, help='Path to the ClinVar file')
    parser.add_argument('--samtools_coverage_file', type=str, help='Path to the samtools coverage file')
    parser.add_argument('--annotation_file', type=str, help='Path to the annotation file')
    parser.add_argument('--output_prefix', type=str, default='output', help='Prefix for output files (default: output)')
    return parser.parse_args()


def combine_files(clinvar_table, coverage_table):
    combined_data = []

    for _, entry in coverage_table.iterrows():
        # Initialize ClinVar-related columns with None
        combined_entry = entry.copy()
        combined_entry['ClinVar_Allele_IDs'] = None
        combined_entry['ClinVar_Variants'] = None
        combined_entry['ClinVar_Review_Status'] = None
        combined_entry['ClinVar_Interpretation'] = None

        matching_rows = clinvar_table[
            (clinvar_table['chrom'] == entry['chrom']) &
            (clinvar_table['position'].between(entry['start'], entry['end']))
        ]

        if not matching_rows.empty:
            for _, matching_row in matching_rows.iterrows():
                combined_entry.update({
                    'ClinVar_Allele_IDs': matching_row['ClinVar_Allele_IDs'],
                    'ClinVar_Variants': matching_row['ClinVar_Variants'],
                    'ClinVar_Review_Status': matching_row['ClinVar_Review_Status'],
                    'ClinVar_Interpretation': matching_row['ClinVar_Interpretation']
                })
                combined_data.append(combined_entry.copy()) 
        else:
            combined_data.append(combined_entry)

    return combined_data



def sort_chromosomes(df):
    def chrom_order(chrom):
        chrom = chrom.replace('chr', '')
        return 23 if chrom == 'X' else (float('inf') if chrom == 'Y' else int(chrom))

    df['chrom_order'] = df['chrom'].apply(chrom_order)
    return df.sort_values(by=['chrom_order', 'start']).drop(columns=['chrom_order'])


def main():
    args = parse_args()

    # Load files
    clinvar_table = pd.read_csv(args.clinvar_file, sep='\t', header=None,
                                names=['chrom', 'position', 'ClinVar_Allele_IDs', 'ClinVar_Variants', 'ClinVar_Review_Status', 'ClinVar_Interpretation'])

    coverage_table = pd.read_csv(args.samtools_coverage_file, sep=',', header=0)
    annotation_table = pd.read_csv(args.annotation_file, sep='\t', header=0)

    # Ensure 'chr' prefix in ClinVar data
    clinvar_table['chrom'] = clinvar_table['chrom'].apply(lambda x: f'chr{x}')


    # Combine ClinVar and coverage data
    clinvar_coverage_df = pd.DataFrame(combine_files(clinvar_table, coverage_table))
    annotation_clinvar_coverage_df = pd.merge(annotation_table, clinvar_coverage_df, on = [ "chrom", "start","end" ], how = "left")
    annotation_clinvar_coverage_df.to_csv("aggregation_test.txt", index = False, sep = "\t")
    # Group and aggregate data
    aggregation = {
        'count_meandepth_gt_20': 'first',
        'count_meanbaseq_gt_20': 'first',
        'count_meanmapq_gt_20': 'first',
        '%_of_samples_with_low_depth': 'first',
        '%_of_samples_with_low_baseq': 'first',
        '%_of_samples_with_low_mapq': 'first',
        'failing_reason': 'first',
        'ClinVar_Allele_IDs': lambda x: ', '.join(sorted(set([str(int(i)) for i in x.dropna() if pd.notnull(i)]))),
        'ClinVar_Variants': lambda x: ', '.join(set(x.dropna())),
        'ClinVar_Review_Status': lambda x: ', '.join(set(x.dropna())),
        'ClinVar_Interpretation': lambda x: ', '.join(set(x.dropna()))
    }

    annotation_clinvar_coverage_df = annotation_clinvar_coverage_df.groupby(['chrom', 'start', 'end']).agg(aggregation).reset_index()

    # Sort chromosomes
    combined_df = sort_chromosomes(annotation_clinvar_coverage_df)

    # Fill missing values and count clinical variants
    combined_df.replace({'': pd.NA, '0': pd.NA}, inplace=True)
    combined_df.fillna('', inplace=True)
    combined_df['#_P/LP_clinical_variants'] = combined_df['ClinVar_Allele_IDs'].apply(
        lambda x: len(x.split(', ')) if isinstance(x, str) and x.strip() else 0)

    # Merge annotation data and write output
    integrated_annotation = pd.merge(annotation_table, combined_df, on=['chrom', 'start', 'end'])
    output_file = f"{args.output_prefix}.integrated_annotation.txt"
    integrated_annotation.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    main()