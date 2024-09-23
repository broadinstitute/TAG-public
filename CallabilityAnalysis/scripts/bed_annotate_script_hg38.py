import argparse
import logging
import pandas as pd
import pybedtools
import warnings
import numpy as np


def parse_args():
    description = "Annotate BED file with Gene/Exon/Transcript information"
    parser = argparse.ArgumentParser(description=description)
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
                        help='Prefix for output files (default: output)', required=False)
    args = parser.parse_args()

    if args.annotation == None:
        raise FileNotFoundError("Annotation DATABASE file not found!")
    if args.bed == None:
        raise FileNotFoundError("Input BED file not found!")
    if args.exome_bed == None:
        raise FileNotFoundError("Input exome BED file not found!")
    if args.gene_bed == None:
        raise FileNotFoundError("Input gene BED file not found!")
    if args.gene_list == None:
        logging.info("Annotating BED file without gene list restrictions")
    if args.output_prefix == None:
        logging.info("Output prefix not specified, using 'output' by default")
    return args


def process_non_coding_bed(unformatted_bed):
    # If the input is empty, return a DataFrame with just the column names
    if unformatted_bed.empty:
        existing_columns = ['chrom', 'start', 'end']
        new_columns = ['interval_size','gene_name', 'exon_number', 'exon_id',
                       'ENSEMBL_transcript_id', 'RefSeq_transcript_id']
        return pd.DataFrame(columns=existing_columns + new_columns)
    
    bed_data = unformatted_bed.copy()
    
    bed_data = bed_data.drop(bed_data.columns[-2:], axis=1)

    existing_columns = ['chrom', 'start', 'end']
    bed_data = bed_data.rename(columns=dict(
        zip(bed_data.columns[:3], existing_columns)))

    bed_data['interval_size'] = bed_data['end'] - bed_data['start'] + 1

    new_columns = ['gene_name', 'exon_number', 'exon_id',
                   'ENSEMBL_transcript_id', 'RefSeq_transcript_id']
    values = (["intergenic"] + [""] * (len(new_columns) - 1))
    bed_data = bed_data.assign(
        **{col: val for col, val in zip(new_columns, values)})

    return bed_data[existing_columns + ['interval_size'] + new_columns]


def process_coding_bed(unformatted_coding_bed, gene_list=None):
    unformatted_coding_bed = unformatted_coding_bed.to_dataframe()
    formatted_coding_bed = pd.DataFrame(columns=['chrom', 'start', 'end', 'interval_size',
                                                 'gene_name', 'exon_number', 'exon_id', 'ENSEMBL_transcript_id', 'RefSeq_transcript_id'])

    for index, row in unformatted_coding_bed.iterrows():
        gene_name, exon_number, exon_id, ENSEMBL_transcript_id, RefSeq_transcript_id = "", "", "", "", ""
        fields = row.values
        chrom, start, end= fields[0], fields[1], fields[2]
        sub_fields = fields[13].split(';')
        for sub_field in sub_fields:
            sub_field = sub_field.strip()
            if sub_field.startswith("gene_name"):
                gene_name = sub_field.split('"')[1]
            if sub_field.startswith("exon_number"):
                exon_number = sub_field.split(" ")[1]
            if sub_field.startswith("exon_id"):
                exon_id = sub_field.split('"')[1]
            if sub_field.startswith("transcript_id"):
                ENSEMBL_transcript_id = sub_field.split('"')[1]
            if sub_field.startswith("db_xref"):
                RefSeq_transcript_id = sub_field.split('"')[1].split(":")[1]
        

        if gene_list is not None and gene_name not in gene_list:
           continue

        interval_size = end - start + 1
        formatted_coding_bed.loc[index] = [chrom, start, end, interval_size, gene_name, exon_number, exon_id, ENSEMBL_transcript_id, RefSeq_transcript_id]

        # Group by chrom, start, end, and concatenate the remaining columns
        grouped_bed = formatted_coding_bed.groupby(['chrom', 'start', 'end', 'interval_size'], as_index=False).agg({
            'gene_name': lambda x: ','.join(sorted(set([str(i) for i in x if i]))),
            'exon_number': lambda x: ','.join(map(str, sorted(set(int(i) for i in x if i)))),
            'exon_id': lambda x: ','.join(sorted(set([str(i) for i in x if i]))),
            'ENSEMBL_transcript_id': lambda x: ','.join(sorted(set([str(i) for i in x if i]))),
            'RefSeq_transcript_id': lambda x: ','.join(sorted(set([str(i) for i in x if i])))
        })

    return grouped_bed

def sort_chromosomes(df):
    def chrom_order(chrom):
        chrom = chrom.replace('chr', '') 
        if chrom in ['X', 'Y']:
            return float('inf') if chrom == 'Y' else 23 
        try:
            return int(chrom) 
        except ValueError:
            return float('inf') 
    
    # Sort by chrom and start
    df['chrom_order'] = df['chrom'].apply(chrom_order)
    sorted_df = df.sort_values(by=['chrom_order', 'start']).drop(columns=['chrom_order'])
    
    return sorted_df


def exonic_bases_calc(input_file, non_exonic_bed):

    # Iterate over the DataFrame rows
    for _, row in non_exonic_bed.iterrows():
        chrom = row['chrom']
        start = int(row['start'])
        end = int(row['end'])
        
        # Calculate the non-exonic interval size
        non_exonic_interval_size = start - end + 1
  
    # Check if non_exonic.bed is empty
    if non_exonic_bed.empty:
        input_file['non_exonic_interval_size'] = 0
    else:
        for idx, row in input_file.iterrows():
            matching_row = [field for field in non_exonic_bed if
                            (field[0] == row['chrom']) and
                            (int(field[1]) >= row['start']) and
                            (int(field[2]) <= row['end'])]
            if matching_row:
                non_exonic_interval_size = int(matching_row[0][5])
                input_file.at[idx,
                                        'non_exonic_interval_size'] = non_exonic_interval_size
            else:
                input_file.at[idx,
                                        'non_exonic_interval_size'] = 0  

    input_file['exonic_interval_size'] = input_file.apply(lambda row: 0 if row['gene_name'] == 'intergenic' else int(row['interval_size'] - row['non_exonic_interval_size']), axis=1)

    input_file.drop(columns=['non_exonic_interval_size'], inplace=True)

    return input_file




def annotate_exon_size(exome_bed, annotated_bed):
    # Convert BedTool to pandas DataFrame
    exome_bed_df = exome_bed.to_dataframe(names=['chrom', 'start', 'end', 'info'])
    exome_bed_df[['gene_name', 'biotype', 'strand', 'region_type']] = exome_bed_df['info'].str.split(':', expand=True)

    exome_bed_df['exonic_interval_size'] = exome_bed_df['end'] - exome_bed_df['start'] + 1

    # Group by gene_name to calculate total exonic bases per gene
    total_exonic_bases_per_gene = exome_bed_df.groupby('gene_name')['exonic_interval_size'].sum().reset_index()
    total_exonic_bases_per_gene['exonic_interval_size'] = total_exonic_bases_per_gene['exonic_interval_size'].astype(int)
    total_exonic_bases_per_gene.columns = ['gene_name', 'total_exonic_bases']

    exome_size_annotated_bed = pd.merge(annotated_bed, total_exonic_bases_per_gene, on='gene_name', how  = "left")
    exome_size_annotated_bed['total_exonic_bases'] = exome_size_annotated_bed['total_exonic_bases'].fillna(0)

    exome_size_annotated_bed['%_exon_impacted'] = exome_size_annotated_bed.apply(
        lambda row: pd.NA if row['total_exonic_bases'] == 0 else round(row['exonic_interval_size'] / row['total_exonic_bases'], 2),
        axis=1
    )
    return exome_size_annotated_bed


def generate_gene_summary(sorted_concatenated_annotated_BED, gene_list=None):
    filtered_bed = sorted_concatenated_annotated_BED[sorted_concatenated_annotated_BED['gene_name'] != 'intergenic']

    # If gene_list_file is provided, filter by genes in that list
    if gene_list:
        filtered_bed = filtered_bed[filtered_bed['gene_name'].isin(gene_list)]
    
    gene_summary = filtered_bed.groupby(['gene_name'], as_index=False).agg({
        'exon_number': lambda x: ','.join(
            map(str, sorted(set(
                int(num) for exon in x for num in str(exon).split(',') if num.strip().isdigit()
            )))
        ),
        'exon_id': lambda x: ','.join(sorted(set([str(i) for i in x if i]))),
        'ENSEMBL_transcript_id': lambda x: ','.join(sorted(set([str(i) for i in x if i]))),
        'RefSeq_transcript_id': lambda x: ','.join(sorted(set([str(i) for i in x if i])))
    })
    
    print(gene_summary)
    return gene_summary

def calculate_coding_and_non_coding_bases(input_bed, total_bases):
    if not input_bed.empty:
        non_coding_base = input_bed[input_bed['gene_name'] == 'intergenic']['interval_size'].sum()
        coding_base = total_bases - non_coding_base
    else:
        # If no intergenic regions, set non-coding bases to 0 and coding bases to total bases
        non_coding_base = 0
        coding_base = total_bases

    return int(non_coding_base), int(coding_base)


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    warnings.filterwarnings('ignore', category=UserWarning, module='pybedtools')
    args = parse_args()

    # Load BED and annotation files
    annotation_file = pybedtools.BedTool(args.annotation)
    bed_file = pybedtools.BedTool(args.bed)
    exome_bed = pybedtools.BedTool(args.exome_bed)

    # Calculate total and coding bases
    total_bases = bed_file.total_coverage()
    coding_bed_file = bed_file.intersect(annotation_file, wa=False, wb=True)
    non_exonic_bed = bed_file.intersect(args.exome_bed, v=True).to_dataframe()
    # intersect to get Coding bases found in MANE gtf 
    coding_bed_to_group_file = bed_file.intersect(annotation_file, wa=True, wb=True)
    # intersect -v to get non-coding bases
    non_coding_bed_to_group = bed_file.intersect(annotation_file, v=True, wa=True).to_dataframe()
    # Process non-coding regions if available
    formatted_non_coding_bed_to_group = process_non_coding_bed(non_coding_bed_to_group)
    # Process coding regions
    gene_list = set(open(args.gene_list).read().splitlines()) if args.gene_list else None
    formatted_coding_bed_to_group = process_coding_bed(coding_bed_to_group_file, gene_list)

    concatenated_annotated_BED = pd.concat([formatted_coding_bed_to_group, formatted_non_coding_bed_to_group])
    sorted_concatenated_annotated_BED = sort_chromosomes(concatenated_annotated_BED) 

    # Calculate Exonic and non-Exonic interval size

    sorted_concatenated_annotated_BED = exonic_bases_calc(sorted_concatenated_annotated_BED, non_exonic_bed)
    fully_annotated_bed = annotate_exon_size(exome_bed,sorted_concatenated_annotated_BED) 
    fully_annotated_bed.to_csv(f"{args.output_prefix}.grouped_by_interval.annotated.txt", sep='\t', index=False)

    # Generate Gene Level Summary
    gene_summary = generate_gene_summary(sorted_concatenated_annotated_BED, gene_list)
    gene_summary.to_csv(f"{args.output_prefix}.grouped_by_gene.txt", sep='\t', index=False)

    genes_involved = len(gene_summary['gene_name'])
    non_coding_base, coding_base = calculate_coding_and_non_coding_bases(fully_annotated_bed, total_bases)

    logging.info(f"This bed file has {non_coding_base} intergenic bases")
    logging.info(f"This bed file has {coding_base} coding bases")
    logging.info(f"This bed file has {genes_involved} genes involved")


    
    
    


    

