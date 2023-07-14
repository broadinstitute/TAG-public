import pybedtools
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Process input files.')
parser.add_argument('--gene_bed', type=str,
                    help='The bed file contains all gene regions with matching gene build with bed file')
parser.add_argument('--bed', type=str, help='Path to the BED file')
parser.add_argument('--grouped_by_gene', type=str,
                    help='Path to previously generated grouped_by_gene annotation file')
parser.add_argument('--gene_list', type=str,
                    help='Path to the gene list file', required=False)
parser.add_argument('--output_prefix', type=str, default='output',
                    help='Prefix for output gene base count file (default: output)')
args = parser.parse_args()


reference_file = pybedtools.BedTool(args.gene_bed)
intersected_regions = pybedtools.BedTool(
    args.bed).intersect(reference_file, wb=True)
total_base_count = pybedtools.BedTool(args.bed).total_coverage()
intersected_regions = intersected_regions.to_dataframe()

intersected_regions = intersected_regions.drop(
    intersected_regions.columns[[3, 4, 5]], axis=1)

intersected_regions.to_csv('gene_file.bed', sep='\t',
                           header=False, index=False)

convert_back_to_bed = pybedtools.BedTool('gene_file.bed').sort()
distinct_genes = {}
total_gene_base_sum = 0
results = []


def generate_gene_base_count(bed, gene_list=None):
    distinct_genes = {}
    results = []
    total_gene_base_sum = 0
    for feature in bed:
        gene_name = feature[3].split(":")[0]
        if args.gene_list:
            with open(args.gene_list, 'r', encoding='UTF-8') as file:
                gene_list = set(file.read().splitlines())
                if gene_name not in gene_list:
                    continue
        if gene_name not in distinct_genes:
            distinct_genes[gene_name] = True
            subset_bed = bed.filter(
                lambda f: f[3].split(":")[0] == gene_name)
            base_count = int(subset_bed.total_coverage())
            results.append((gene_name, base_count))
            total_gene_base_sum += base_count
    results = sorted(results, key=lambda x: x[0])  # Sort by gene name
    gene_base_count = pd.DataFrame(
        results, columns=["gene_name", "base_count"])
    return gene_base_count


gene_base_count = generate_gene_base_count(convert_back_to_bed)

gene_name_list = gene_base_count['gene_name'].tolist()
num_elements = int(len(gene_name_list)) - 1
print(num_elements)


def generate_total_gene_base_count(gene_bed, gene_list=None):
    distinct_genes = {}
    results = []
    gene_bed = pd.read_csv(gene_bed, sep="\t")
    all_genes = gene_bed.iloc[:, 3].str.split(':').str[0]
    gene_bed_genes = set(all_genes.unique())
    total_gene_base_count = pd.DataFrame(
        columns=['gene_name', 'total_base_count'])
    for gene_name in gene_name_list:
        if gene_name not in gene_bed_genes:
            continue
        if gene_name not in distinct_genes:
            distinct_genes[gene_name] = True
            gene_bed_df = pybedtools.BedTool(args.gene_bed)
            subset_bed = gene_bed_df.filter(
                lambda f: f[3].split(":")[0] == gene_name)
            total_base_count = int(subset_bed.total_coverage())
            results.append((gene_name, total_base_count))
    for gene_name, total_base_count in results:
        total_gene_base_count = pd.concat([
            total_gene_base_count,
            pd.DataFrame(
                {'gene_name': [gene_name], 'total_base_count': [total_base_count]})
        ], ignore_index=True)
    return total_gene_base_count


total_gene_base_count = generate_total_gene_base_count(args.gene_bed)


merged_base_info = pd.merge(
    gene_base_count, total_gene_base_count, on="gene_name", how="left")

grouped_by_gene = pd.read_csv(args.grouped_by_gene, sep='\t', header=0)

merged_gene_info = pd.merge(
    grouped_by_gene, merged_base_info, on="gene_name", how="inner")

merged_gene_info['percentage'] = merged_gene_info['base_count'] / \
    merged_gene_info['total_base_count']
header = ['gene_name', 'exon_number', 'exon_id',
          'transcript_id', 'base_count', 'total_base_count', 'percentage']
merged_gene_info.columns = header

merged_gene_info.to_csv(
    args.output_prefix+'.gene_base_count.txt', index=False, sep='\t')
