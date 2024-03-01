# Running this script is a prerequisite for demux_CODEC wdl.
# This script is taking an Excel file provided by collaborator
# with Barcode assigned to each sample and a default Index reference file,
# to generate a sample_sheet_L00{lane}.csv file for each lane.


import sys
import pandas as pd


if len(sys.argv) != 3:
    print("Usage: python3 prep_codec_metadata.py <Metadata Excel file> <Index CSV file>")
    sys.exit(1)

metadata_excel_file = sys.argv[1]
index_csv_file = sys.argv[2]

xlsx_df = pd.read_excel(metadata_excel_file)
index_df = pd.read_csv(index_csv_file)

merged_df = pd.merge(
    xlsx_df, index_df, left_on='CODEC index', right_on='index')

for lane, group in merged_df.groupby('lanes'):

    output_df = group[['submission_id', 'IndexBarcode1', 'IndexBarcode2']]
    output_df.columns = ['SampleName', 'IndexBarcode1', 'IndexBarcode2']

    output_df.to_csv(f'sample_sheet_L00{lane}.csv', index=False)
