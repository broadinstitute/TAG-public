import pandas as pd
import sys


file_paths = sys.argv[1].split()

if len(sys.argv) > 2:
    threshold = int(sys.argv[2])
else:
    print("Threshold not provided. Using default value of 80.")
    threshold = 80

counts_dict = {}
columns_to_count = ['meandepth', 'meanbaseq', 'meanmapq']
dfs = []
for file_path in file_paths:
    df = pd.read_csv(file_path, sep='\t', usecols=[
        '#rname', 'startpos', 'endpos', 'meandepth', 'meanbaseq', 'meanmapq'],
        comment=None, header=0)
    dfs.append(df)
# Initialize counts_dict with all keys and zero counts
for df in dfs:
    for index, row in df.iterrows():
        key = (row['#rname'], row['startpos'], row['endpos'])
        if key not in counts_dict:
            counts_dict[key] = {col: 0 for col in columns_to_count}

# Iterate over each line and each file to count
for df in dfs:
    for index, row in df.iterrows():
        rname = row['#rname']
        startpos = row['startpos']
        endpos = row['endpos']
        for column in columns_to_count:
            column_gt_20 = row[column] > 20
            if column_gt_20:
                key = (rname, startpos, endpos)
                counts_dict[key][column] += 1

results_df = pd.DataFrame.from_dict(counts_dict, orient='index').reset_index()

results_df.rename(columns={
    'level_0': 'chrom',
    'level_1': 'start',
    'level_2': 'end',
    'meandepth': 'count_meandepth_gt_20',
    'meanbaseq': 'count_meanbaseq_gt_20',
    'meanmapq': 'count_meanmapq_gt_20'
}, inplace=True)


results_df['chrom'] = results_df['chrom'].str.replace('chr', '')
# Initial conversion
results_df['chrom'] = results_df['chrom'].apply(
    lambda x: '999' if x == 'X' else ('1000' if x == 'Y' else x))

# Sort the DataFrame by genomic coordinates
results_df = results_df.sort_values(by=['chrom', 'start'])

results_df['chrom'] = results_df['chrom'].apply(
    lambda x: 'chrX' if x == '999' else ('chrY' if x == '1000' else f'chr{x}'))

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
               >= threshold, 'failing_reason'] += 'depth,'
results_df.loc[results_df['%_of_samples_with_low_baseq']
               >= threshold, 'failing_reason'] += 'baseq,'
results_df.loc[results_df['%_of_samples_with_low_mapq']
               >= threshold, 'failing_reason'] += 'mapq,'
results_df['failing_reason'] = results_df['failing_reason'].str.rstrip(',')
results_df.to_csv('samtools_coverage_summary.txt', index=False)