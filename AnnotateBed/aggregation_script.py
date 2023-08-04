import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='Add input files.')
parser.add_argument('--clinvar_file', type=str,
                    help='Path to the clinvar file')
parser.add_argument('--samtools_coverage_file', type=str,
                    help='Path to the samtools coverage file')
parser.add_argument('--annotation_file', type=str,
                    help='Path to the annotation file')
parser.add_argument('--output_prefix', type=str, default='output',
                    help='Prefix for output files (default: output)')
args = parser.parse_args()


def combine_files(clinvar_table, coverage_table):
    combined_data = []
    for _, entry in coverage_table.iterrows():
        chrom, start, end, count_meandepth_gt_20, count_meanbaseq_gt_20, count_meanmapq_gt_20, pct_of_samples_with_low_depth, pct_of_samples_with_low_baseq, pct_of_samples_with_low_mapq, failing_reason = entry['chrom'], entry[
            'start'], entry['end'], entry['count_meandepth_gt_20'], entry['count_meanbaseq_gt_20'], entry['count_meanmapq_gt_20'], entry['%_of_samples_with_low_depth'], entry['%_of_samples_with_low_baseq'], entry['%_of_samples_with_low_mapq'], entry['failing_reason']
        matching_rows = clinvar_table[
            (clinvar_table['chrom'] == chrom) &
            (clinvar_table['position'] >= int(start)) &
            (clinvar_table['position'] <= int(end))
        ]
        if not matching_rows.empty:
            if not matching_rows.empty:
                for _, matching_row in matching_rows.iterrows():
                    combined_entry = entry.copy()
                    combined_entry['ClinVar_Allele_IDs'] = matching_row['ClinVar_Allele_IDs']
                    combined_entry['ClinVar_Variants'] = matching_row['ClinVar_Variants']
                    combined_entry['ClinVar_Review_Status'] = matching_row['ClinVar_Review_Status']
                    combined_entry['ClinVar_Interpretation'] = matching_row['ClinVar_Interpretation']
                    combined_data.append(combined_entry)
        else:
            combined_data.append(entry)
    return combined_data


aggregation = {
    'count_meandepth_gt_20': 'first',
    'count_meanbaseq_gt_20': 'first',
    'count_meanmapq_gt_20': 'first',
    '%_of_samples_with_low_depth': 'first',
    '%_of_samples_with_low_baseq': 'first',
    '%_of_samples_with_low_mapq': 'first',
    'failing_reason': 'first',
    'ClinVar_Allele_IDs': lambda x: ', '.join(((x.dropna()))),
    'ClinVar_Variants': lambda x: ', '.join(set(x.dropna())),
    'ClinVar_Review_Status': lambda x: ', '.join(set(x.dropna())),
    'ClinVar_Interpretation': lambda x: ', '.join(set(x.dropna()))
}


def main():
    clinvar_table = pd.read_csv(
        args.clinvar_file, sep='\t', header=None, names=['chrom', 'position', 'ClinVar_Allele_IDs', 'ClinVar_Variants', 'ClinVar_Review_Status', 'ClinVar_Interpretation']
    )
    coverage_table = pd.read_csv(
        args.samtools_coverage_file, sep='\t', header=0)
    annotation_table = pd.read_csv(
        args.annotation_file, sep='\t', header=0)
    clinvar_table['chrom'] = clinvar_table['chrom'].apply(lambda x: f'chr{x}')
    combined_data = combine_files(clinvar_table, coverage_table)
    combined_data_df = pd.DataFrame(combined_data, columns=coverage_table.columns.tolist() +
                                    clinvar_table.columns[2:].tolist())
    print(combined_data_df.iloc[:, -4])
    combined_data_df['ClinVar_Allele_IDs'] = pd.to_numeric(combined_data_df['ClinVar_Allele_IDs'],
                                                           errors='coerce').fillna(0).astype(int).astype(str)
    combined_data_df = combined_data_df.groupby(
        ['chrom', 'start', 'end']).agg(aggregation).reset_index()
    combined_data_df['chrom'] = combined_data_df['chrom'].str.replace(
        r'^chr', '', regex=True).astype(int)
    combined_data_df['#_P/LP_clinical_variants'] = combined_data_df['ClinVar_Allele_IDs'].apply(
        lambda x: len(x.split(', ')) if x != 'NaN' and isinstance(
            x, str) and x.strip() else 0
    )
    # Sort the DataFrame by genomic coordinates (first column, then second column)
    combined_data_df = combined_data_df.sort_values(by=['chrom', 'start'])
    combined_data_df['chrom'] = combined_data_df['chrom'].apply(
        lambda x: f'chr{x}')
    combined_data_df.replace({'': pd.NA, '0': pd.NA},
                             inplace=True)
    combined_data_df.fillna('', inplace=True)
    combined_data_df['#_P/LP_clinical_variants'] = combined_data_df['ClinVar_Allele_IDs'].apply(
        lambda x: len(x.split(', ')) if x != 'NaN' and isinstance(
            x, str) and x.strip() else 0
    )
    print(combined_data_df)
    merged_df = pd.concat(
        [annotation_table, combined_data_df.iloc[:, 3:]], axis=1)
    output_file_name = args.output_prefix + '.integrated_annotation.txt'
    merged_df.to_csv(output_file_name, sep='\t', index=False)


if __name__ == "__main__":
    main()
