import pandas as pd
import sys


file_paths = sys.argv[1].split()

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
