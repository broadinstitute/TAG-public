"""
AnnotateBed

This script provides annotation for a given bed file, annotation GTF and a bed file for exome with matching gene build.

Author: Stella Li
Date: 2023-05-30
"""

import argparse
import pandas as pd
import pybedtools

# Create an argument parser
parser = argparse.ArgumentParser(description='Process input files.')

# Add the command-line arguments
parser.add_argument('--annotation', type=str,
                    help='Path to the annotation file')
parser.add_argument('--bed', type=str, help='Path to the BED file')
parser.add_argument('--gene_bed', type=str,
                    help='The bed file contains all gene regions with matching gene build with bed file')
parser.add_argument('--exome_bed', type=str,
                    help='The bed file contains all exon regions with matching gene build with bed file')
parser.add_argument('--gene_list', type=str,
                    help='Path to the gene list file', required=False)
parser.add_argument('--output_prefix', type=str, default='output',
                    help='Prefix for output files (default: output)')
args = parser.parse_args()

# Use the provided file paths
annotation_file = pybedtools.BedTool(args.annotation)
bed_file = pybedtools.BedTool(args.bed)
# count bases in given bed file
total_bases = bed_file.total_coverage()

# intersect to get the bases in GTF. We want wa=False in order to count bases of genes and we want wa=True to get grouped annotation
coding_bed_file = bed_file.intersect(annotation_file, wa=False, wb=True)
coding_bed_to_group_file = bed_file.intersect(
    annotation_file, wa=True, wb=True)
# intersect to get the bases not in GTF, that is to say, non-coding bases
non_coding_bed = bed_file.intersect(annotation_file, v=True).to_dataframe()
non_coding_bed_to_group = bed_file.intersect(
    annotation_file, v=True, wa=True).to_dataframe()


def process_non_coding_bed_to_group_file(unformatted_bed):
    bed_data = unformatted_bed.copy()
    bed_data = bed_data.drop(bed_data.columns[-2:], axis=1)
    existing_columns = ['chrom', 'start', 'end']
    bed_data = bed_data.rename(columns=dict(
        zip(bed_data.columns[:3], existing_columns)))
    new_columns = ['gene_name', 'exon_number', 'transcript_id']
    values = [""] * (len(new_columns) - 1)
    bed_data = bed_data.assign(
        **{col: val for col, val in zip(new_columns, values)})
    return bed_data


def process_non_coding_bed_file(unformatted_bed):
    bed_data = unformatted_bed.copy()
    bed_data = bed_data.drop(bed_data.columns[-2:], axis=1)
    existing_columns = ['chrom', 'start', 'end']
    bed_data = bed_data.rename(columns=dict(
        zip(bed_data.columns[:3], existing_columns)))
    new_columns = ['function', 'gene_type', 'gene_name',
                   'exon_number', 'exon_id', 'transcript_id']
    values = ["intergenic"] + [""] * (len(new_columns) - 1)
    bed_data = bed_data.assign(
        **{col: val for col, val in zip(new_columns, values)})
    return bed_data


# Process non-coding bed file (prepared for ungrouped output)
formatted_non_coding_bed = pd.DataFrame()
if not non_coding_bed.empty:
    formatted_non_coding_bed = process_non_coding_bed_file(non_coding_bed)
    formatted_non_coding_bed_to_group = process_non_coding_bed_to_group_file(
        non_coding_bed_to_group)


def process_coding_bed_file(unformatted_coding_bed, gene_list=None):
    unformatted_non_coding_bed = unformatted_coding_bed.to_dataframe()
    formatted_coding_bed = pd.DataFrame(columns=['chrom', 'start', 'end', 'function',
                                                 'gene_type', 'gene_name', 'exon_number', 'exon_id', 'transcript_id'])
    for index, row in unformatted_non_coding_bed.iterrows():
        gene_name, gene_type, exon_number, exon_id, transcript_id = "", "", "", "", ""
        fields = row.values
        chrom, start, end, function = fields[5], fields[8], fields[9], fields[7]
        sub_fields = fields[13].split(';')
        for sub_field in sub_fields:
            sub_field = sub_field.strip()
            if sub_field.startswith("gene_type"):
                gene_type = sub_field.split('"')[1]
            if sub_field.startswith("gene_name"):
                gene_name = sub_field.split('"')[1]
            if sub_field.startswith("exon_number"):
                exon_number = sub_field.split(" ")[1]
            if sub_field.startswith("exon_id"):
                exon_id = sub_field.split('"')[1]
            if sub_field.startswith("transcript_id"):
                transcript_id = sub_field.split('"')[1]
        if gene_list is not None and gene_name not in gene_list:
            continue
        formatted_coding_bed.loc[index] = [chrom, start, end, function,
                                           gene_type, gene_name, exon_number, exon_id, transcript_id]
    return formatted_coding_bed


gene_list = None
if args.gene_list:
    with open(args.gene_list, 'r', encoding='UTF-8') as file:
        gene_list = set(file.read().splitlines())


formatted_coding_bed = process_coding_bed_file(coding_bed_file, gene_list)

# concatenate coding and non-coding files into a single annotation file if non-coding file is not empty
if not non_coding_bed.empty:
    # concatenate non coding regions and coding regions
    combined_file = pd.concat(
        [formatted_non_coding_bed, formatted_coding_bed])
    combined_file.to_csv('ungrouped_annotation.txt', sep='\t',
                         header=False, index=False)
else:
    # take only infos from coding file
    header = ['chrom', 'start', 'end', 'function', 'gene_type',
              'gene_name', 'exon_number', 'exon_id', 'transcript_id']
    formatted_coding_bed.to_csv('ungrouped_annotation.txt',
                                sep='\t', header=header, index=False)

# Sort by genomic coordinates and remove annotations for 'gene'


def sort_filter_annotation_file(input_file, output_file):
    with open(input_file, 'r', encoding='utf-8') as input_file:
        lines = input_file.readlines()
        header = lines[0]
        data = lines[1:]
        sorted_data = sorted(data, key=lambda x: (
            int(x.split('\t')[0][3:]) if x.split('\t')[
                0][3:].isdigit() else float('inf'),
            int(x.split('\t')[1])))
        seen_lines = set()
        filtered_lines = []
        for line in sorted_data:
            if line not in seen_lines and line.split('\t')[3] != 'gene':
                filtered_lines.append(line)
                seen_lines.add(line)
        filtered_content = header + ''.join(filtered_lines)
        with open(output_file, 'w', encoding='utf-8') as output_file:
            output_file.write(filtered_content)


sort_filter_annotation_file(
    'ungrouped_annotation.txt', 'tmp_ungrouped_annotation.txt')


def sort_internal_elements(row):
    row['transcript_id'] = ','.join(
        sorted(x for x in row['transcript_id'].split(',')))
    exon_numbers = [int(x) for x in row['exon_number'].split(
        ',') if x.strip()]  # Skip empty strings
    sorted_exon_numbers = sorted(exon_numbers)
    row['exon_number'] = ','.join(str(x) for x in sorted_exon_numbers)
    return row


input_file = pd.read_csv('tmp_ungrouped_annotation.txt', sep='\t', header=0)

input_file.columns = ['chrom', 'start', 'end', 'function', 'gene_type',
                      'gene_name', 'exon_number', 'exon_id', 'transcript_id']
input_file['exon_number'] = pd.to_numeric(
    input_file['exon_number'], errors='coerce').fillna(0)
input_file['exon_number'] = input_file['exon_number'].astype(int)
input_file = input_file[input_file['exon_number'] != 0]
input_file['exon_number'] = input_file['exon_number'].astype(str)


grouped_by_gene = input_file.groupby(['gene_name']).agg({'exon_number': lambda x: ', '.join(
    set(x.dropna())), 'exon_id': lambda x: ', '.join(set(x.dropna())), 'transcript_id': lambda x: ', '.join(set(x.dropna()))}).reset_index()

grouped_by_gene[['exon_number', 'exon_id', 'transcript_id']] = grouped_by_gene[['exon_number',
                                                                                'exon_id', 'transcript_id']].apply(sort_internal_elements, axis=1)

grouped_by_gene_file = args.output_prefix + '.grouped_by_gene.txt'
grouped_by_gene.to_csv(grouped_by_gene_file, sep='\t', index=False)

# write ungrouped annotation file
ungrouped_annotation_file = args.output_prefix + '.ungrouped.annotated.txt'
input_file.to_csv(ungrouped_annotation_file, sep='\t', index=False)

# The code below is to generate one line per interval grouped file


def process_coding_bed_to_group_file(unformatted_coding_bed, gene_list=None):
    unformatted_non_coding_bed = unformatted_coding_bed.to_dataframe()
    formatted_coding_bed = pd.DataFrame(
        columns=['chrom', 'start', 'end', 'gene_name', 'exon_number', 'transcript_id'])
    for index, row in unformatted_non_coding_bed.iterrows():
        gene_name, exon_number, transcript_id = "", "", ""
        fields = row.values
        chrom, start, end, function = fields[0], fields[1], fields[2], fields[7]
        sub_fields = fields[13].split(';')
        for sub_field in sub_fields:
            sub_field = sub_field.strip()
            if sub_field.startswith("gene_name"):
                gene_name = sub_field.split('"')[1]
            if sub_field.startswith("exon_number"):
                exon_number = sub_field.split(" ")[1]
            if sub_field.startswith("transcript_id"):
                transcript_id = sub_field.split('"')[1]
        if gene_list is not None and gene_name not in gene_list:
            continue
        formatted_coding_bed.loc[index] = [chrom, start, end,
                                           gene_name, exon_number, transcript_id]
    return formatted_coding_bed


formatted_coding_bed_to_group = process_coding_bed_to_group_file(
    coding_bed_to_group_file, gene_list)

# concatenate coding and non-coding files into a single annotation file if non-coding file is not empty
if not non_coding_bed.empty:
    # concatenate non coding regions and coding regions
    combined_file = pd.concat(
        [formatted_non_coding_bed_to_group, formatted_coding_bed_to_group])
    combined_file.to_csv('grouped_annotation.txt', sep='\t',
                         header=0, index=False)
else:
    # take only infos from coding file
    header = ['chrom', 'start', 'end', 'gene_name',
              'exon_number', 'transcript_id']
    formatted_coding_bed_to_group.to_csv('grouped_annotation.txt',
                                         sep='\t', header=header, index=False)

grouped_by_gene[['exon_number', 'transcript_id']] = grouped_by_gene[[
    'exon_number', 'transcript_id']].apply(sort_internal_elements, axis=1)

sort_filter_annotation_file(
    'grouped_annotation.txt', 'tmp_grouped_annotation.txt')

input_file = pd.read_csv('tmp_grouped_annotation.txt', sep='\t', header=0)
input_file.columns = ['chrom', 'start', 'end',
                      'gene_name', 'exon_number', 'transcript_id']
input_file['exon_number'] = pd.to_numeric(
    input_file['exon_number'], errors='coerce').fillna(0).astype(int)
input_file['exon_number'] = input_file['exon_number'].astype(str)

grouped_input_file = input_file.groupby(['chrom', 'start', 'end']).agg({'gene_name': lambda x: ', '.join(set((x.dropna()))), 'exon_number': lambda x: ', '.join(
    set(x.dropna())), 'transcript_id': lambda x: ', '.join(set(x.dropna()))}).reset_index()
columns_to_clean = ['exon_number', 'transcript_id']

for column in columns_to_clean:
    grouped_input_file[column] = grouped_input_file[column].str.split(',').apply(lambda x: [value.strip(
    ) for value in x if value.strip() not in ('nan', '<NA>', '0')]).str.join(',')

# sort the value for these columns: ['exon_number', 'exon_id','transcript_id']
grouped_input_file[['exon_number', 'transcript_id']] = grouped_input_file[['exon_number',
                                                                           'transcript_id']].apply(sort_internal_elements, axis=1)

# code below is to generate exon base count


def calculate_interval_size(start, end):
    return end - start


exome_bed_file = pybedtools.BedTool(args.exome_bed)
non_exonic_bed_file = bed_file.intersect(exome_bed_file, v=True)
non_exonic_bed_file.saveas('non_exonic.bed')
count_exon_bases = []
with open(args.exome_bed, 'r') as file:
    for line in file:
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        name = fields[3]
        gene_name = name.split(':')[0]
        sub_interval_size = calculate_interval_size(start, end)
        count_exon_bases.append(
            (chrom, start, end, gene_name, sub_interval_size))


# Calculate the exon size for each unique gene_name
exon_sizes = {}
for entry in count_exon_bases:
    gene_name = entry[3]
    sub_interval_size = entry[4]
    if gene_name not in exon_sizes:
        exon_sizes[gene_name] = sub_interval_size
    else:
        exon_sizes[gene_name] += sub_interval_size


file1_data = []
with open('non_exonic.bed', 'r') as file1:
    for line in file1:
        fields = line.strip().split('\t')
        non_exonic_interval_size = calculate_interval_size(
            int(fields[1]), int(fields[2]))
        file1_data.append((fields[0], int(fields[1]),
                          int(fields[2]), fields[3], fields[4], non_exonic_interval_size))


grouped_input_file.insert(
    3, 'interval_size', grouped_input_file['end'] - grouped_input_file['start'])
for idx, row in grouped_input_file.iterrows():
    matching_row = [field for field in file1_data if
                    (field[0] == row['chrom']) and
                    (int(field[1]) <= row['start']) and
                    (int(field[2]) >= row['end'])]
    if matching_row:
        non_exonic_interval_size = int(matching_row[0][5])
        grouped_input_file.at[idx,
                              'non_exonic_interval_size'] = non_exonic_interval_size
    else:
        grouped_input_file.at[idx, 'non_exonic_interval_size'] = 0

grouped_input_file['size_of_exon'] = grouped_input_file['interval_size'] - \
    grouped_input_file['non_exonic_interval_size']
grouped_input_file['total_exon_size'] = grouped_input_file['gene_name'].map(
    exon_sizes)
grouped_input_file['%_exon_impacted_for'] = grouped_input_file['size_of_exon'] / \
    grouped_input_file['total_exon_size']
grouped_input_file.drop(
    columns=['non_exonic_interval_size', 'total_exon_size'], inplace=True)
output_grouped_file = args.output_prefix + '.grouped_by_interval.annotated.txt'

grouped_input_file.to_csv(output_grouped_file, sep='\t', index=False)

with open(output_grouped_file, 'r', encoding='utf-8') as input_file:
    lines = input_file.readlines()
    header = lines[0]
    data = lines[1:]
    # Sort the data by genomic coordinates
    sorted_data = sorted(data, key=lambda x: (
        int(x.split('\t')[0][3:]) if x.split('\t')[
            0][3:].isdigit() else float('inf'),
        int(x.split('\t')[1])))
    sorted_content = header + ''.join(sorted_data)

with open(output_grouped_file, 'w', encoding='utf-8') as output_file:
    output_file.writelines(sorted_content)

# The code below is to generate a file with gene base count file
if not non_coding_bed.empty:
    non_coding_bed['base_count'] = pd.to_numeric(
        non_coding_bed.iloc[:, 2], errors='coerce') - pd.to_numeric(non_coding_bed.iloc[:, 1], errors='coerce')
    intergenic_base = int(non_coding_bed['base_count'].sum())
    coding_base = total_bases - intergenic_base
else:
    intergenic_base = 0
    coding_base = total_bases

with open("intergenic_base_count.txt", "w") as file:
    file.write(f"{intergenic_base}")
with open("coding_base_count.txt", "w") as file:
    file.write(f"{coding_base}")
