# cnvArrayProber

## Overview

The `cnvArrayProber` is designed to analyze CNV (Copy Number Variation) intervals from a BED file and map probe information from two array support files (CytoSNP-850K and GDA). It generates a comprehensive XLSX file containing the intervals and the number of probes they contain in each array. Additionally, the script outputs a PDF file with detailed plots illustrating the locations of these probes in the CytoSNP and GDA arrays.

This script was developed upon a request from Greg Nakashian in [TAG1994](https://github.com/broadinstitute/TAG/issues/1994).

## Features

1. **Input Processing:**
   - **BED File:** The script reads CNV intervals from a specified BED file.
   - **Array Support Files:** It also processes two array support files, CytoSNP-850K and GDA, to gather probe information.

2. **Data Analysis:**
   - **Interval Analysis:** The [cnvArrayProber](https://dockstore.org/workflows/github.com/broadinstitute/TAG-public/cnvArrayProber:array_prober_yg?tab=info) WDL analyzes the CNV intervals to determine the number of probes from each array (CytoSNP-850K and GDA) that fall within each interval.
   
3. **Output Generation:**
   - **XLSX File:** A xlsx file is generated, containing the CNV intervals and the corresponding count of probes from each array.
   - **PDF File:** A PDF file is produced, featuring plots that visually represent the locations of the probes within the CytoSNP and GDA arrays.

## Usage

To use the `cnvArrayProber` WDL, follow these steps:

1. **Prepare Input Files:**
   - Ensure you have a BED file containing the CNV intervals.
```
chr2	97220584	130400286
chr4	143920938	144022444
chr9	33140790	33261063
  ```
	
You can get CytoSNP-850K and GDA array support files using the following gcloud link. (**Note: Ensure you are using consistent genome build for those input**)


| GDACyto_hg19_SupportFile                                                   | CytoSNP_850k_v1_4_hg38_SupportFile                                                                      | GDACyto_hg38_SupportFile                                                             | CytoSNP_850k_v1_4_hg19_SupportFile                                                                      |
|:-------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------|
| gs://fc-d2c7d48c-9433-4a1f-bdeb-100265b01a63/GDA_SupportFile/GDACyto_20047166_A1.csv | gs://fc-d2c7d48c-9433-4a1f-bdeb-100265b01a63/CytoSNP-850Kv1-4_SupportFile/CytoSNP-850Kv1-4_iScan_B2.csv | gs://fc-d2c7d48c-9433-4a1f-bdeb-100265b01a63/GDA_SupportFile/GDACyto_20047166_A2.csv | gs://fc-d2c7d48c-9433-4a1f-bdeb-100265b01a63/CytoSNP-850Kv1-4_SupportFile/CytoSNP-850Kv1-4_iScan_B1.csv |



2. **Execute the Workflow:**
   - Execute the `cnvArrayProber` WDL with inputs defined by data table.
   - The WDL will process the files and generate the output XLSX and PDF files.

3. **Review Outputs:**

- **XLSX File:** A Exel file with the following information in two separate sheets for CytoSNP-850K and GDA arrays:

CytoSNP850K:

|                          |   left_padding |   interval |   right_padding |
|:-------------------------|---------------:|-----------:|----------------:|
| chr2:97220584-130400286  |            239 |       8051 |             171 |
| chr4:143920938-144022444 |              0 |          4 |               1 |
| chr9:33140790-33261063   |              2 |         38 |               2 |

GDA:

|                          |   left_padding |   interval |   right_padding |
|:-------------------------|---------------:|-----------:|----------------:|
| chr2:97220584-130400286  |            897 |      19284 |             683 |
| chr4:143920938-144022444 |              2 |         41 |               1 |
| chr9:33140790-33261063   |              9 |         81 |              18 |


- **PDF File:** A PDF document with plots showing:
  - The distribution of CytoSNP-850K probes within each interval.
  - The distribution of GDA probes within each interval.


## Development and Contributions


The script was developed by Yueyao Gao (gaoyueya@broadinstitute.org) in response to a request from Greg Nakashian in [TAG1994](https://github.com/broadinstitute/TAG/issues/1994). Contributions and further improvements are welcome. Please refer to the TAG repo for more information.

