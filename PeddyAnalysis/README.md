# Peddy Analysis Pipeline

This PeddyAnalysis workflow is designed to use gVCFs to generate relatedness statistics within familiy (Can process multiple families at a time and should submit at sample set level). The pipeline helps analyze and visualize familial relationships and concordance with known pedigree information using Peddy.

**NOTE: For prediction accuracy, gVCFs are encouraged for this analysis.**

## Inputs
- Array[String] sample_ids: List of sample IDs.
- Array[String] family_ids: List of family IDs corresponding to each sample.
- Array[String] gvcf_paths: Paths to gVCF files for each sample.
- Array[String] gvcf_index_paths: Paths to gVCF index files.
- Array[String] pedigrees: Pedigree information for each sample. [Father/Mother/Proband].
- Array[String] reported_sexes: Reported sexes of each sample. [Female/Male]
- File interval_list: File defining the genomic intervals for evaluation.
- File reference_fasta: Reference genome in FASTA format.
- File reference_fasta_index: Index file for the reference genome.
- File reference_dict: Dictionary file for the reference genome.

## Outputs

- File family_composition_log: Log file summarizing family compositions (trios, duos, single-sample families).
- Int num_single_sample_families: Number of families with only one sample. (Family with one sample only will not be analyzed) 
- Array[String] processed_families: List of processed families.
- Array[String] unprocessed_family_ids: List of families that were not processed.
- Array[File] merged_gvcf_files: The genotyped gVCF files for each family.
- Array[File] merged_gvcf_indices: The indices for the genotyped gVCFs.
- Array[File] family_info_files: Information for family decompostion.
- Array[File] updated_fam_file: Updated .fam files with corrected family info.
- Float concordance: The concordance of pedigree predictions compared to the known family information.
- File merged_peddy_results: Merged Peddy results for all families.
- File all_family_peddy_prediction_plot: Plot of relatedness and pedigree prediction for all families processed.

## Workflow Diagram


![workflow](https://github.com/user-attachments/assets/b963f544-7516-41a6-a924-802300771437)


## Example Delivery


#### 1. Summarized Peddy Statistics


| family_id | sample_a | sample_b | rel       | hets_a | hets_b | shared_hets | ibs0 | ibs2 | n    | pedigree_parents | pedigree_relatedness | predicted_parents | parent_error | sample_duplication_error | rel_difference |
| --------- | -------- | -------- | --------- | ------ | ------ | ----------- | ---- | ---- | ---- | ---------------- | -------------------- | ----------------- | ------------ | ------------------------ | -------------- |
| SF0429362 | SDSM-TX  | SDSM-TY  | 0.4925    | 3624   | 3708   | 1795        | 5    | 3075 | 6822 | TRUE             | 0.5                  | TRUE              | FALSE        | FALSE                    | 0.00745        |
| SF0026276 | SDSM-WN  | SDSM-WO  | 0.507     | 3781   | 3722   | 1897        | 5    | 3169 | 6867 | TRUE             | 0.5                  | TRUE              | FALSE        | FALSE                    | \-0.006985     |
| SF0051349 | SDSM-U6  | SDSM-U7  | \-0.02671 | 2546   | 2656   | 972         | 520  | 2537 | 6303 | FALSE            | 0                    | FALSE             | FALSE        | FALSE                    | 0.02671        |
| SF0051349 | SDSM-U6  | SDSM-U8  | 0.4921    | 2546   | 2647   | 1255        | 1    | 3624 | 6303 | TRUE             | 0.5                  | TRUE              | FALSE        | FALSE                    | 0.007855       |
| SF0051349 | SDSM-U7  | SDSM-U8  | 0.5006    | 2656   | 2647   | 1339        | 7    | 3690 | 6315 | TRUE             | 0.5                  | TRUE              | FALSE        | FALSE                    | \-0.0005667    |
| SF0065096 | SDSM-W4  | SDSM-WB  | \-0.01655 | 3565   | 3590   | 1107        | 583  | 2176 | 7698 | FALSE            | 0                    | FALSE             | FALSE        | FALSE                    | 0.01655        |
| SF0065096 | SDSM-W4  | SDSM-WC  | 0.4884    | 3565   | 3640   | 1743        | 1    | 3979 | 7698 | TRUE             | 0.5                  | TRUE              | FALSE        | FALSE                    | 0.01164        |
| SF0065096 | SDSM-WB  | SDSM-WC  | 0.5145    | 3590   | 3640   | 1861        | 7    | 4196 | 7699 | TRUE             | 0.5                  | TRUE              | FALSE        | FALSE                    | \-0.01448      |
| SF0417467 | SDSM-TI  | SDSM-TL  | 0.5146    | 3739   | 3725   | 1931        | 7    | 4258 | 7859 | TRUE             | 0.5                  | TRUE              | FALSE        | FALSE                    | \-0.01463      |
| SF0417467 | SDSM-TI  | SDSM-VS  | 0.01341   | 3739   | 3653   | 1141        | 546  | 2195 | 7859 | FALSE            | 0                    | FALSE             | FALSE        | FALSE                    | \-0.01341      |


#### 2. Concordance


![concordance-check](https://github.com/user-attachments/assets/55e66e6d-da31-4860-8c03-6928ad4a7a03)


