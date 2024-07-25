# CramToBam Workflow
#### CramToBam
This task converts the inputted Cram file into a Bam file. It requires the input Cram file along with its reference genome dict, fasta file, and the fasta.fai reference index file. The outputs include the Bam file and the Bam_index_file (.bam.bai).
#### ValidateBam
The second task simply validates that the Bam File has been correctly created, and it outputs a summary .txt file with any warnings or errors that it found in the bam file.


## Usage with Terra

1. This CramToBam Conversion task is currently only available via Terra. You can find the WDL [here](https://app.terra.bio/#workspaces/broadtagteam/tag2010_PECGS_cvg_drop/workflows/TAG-PUBLIC-BUGFIXED-WDLS/FixedCramToBam). 


2. We will now input some test data to exhibit the functionality of the WDL. For the test inputs, we have the .cram file for Sample (INPUT SAMPLE HERE). And since this sample is referenced against the hg19 reference genome, we will use the hg19 reference files.

3. To input the aforementioned files into the WDL, paste the following strings into the respective input variables:

###### Sample Data:
```
"SAMPLE .cram FILE INPUT HERE"
```
###### Reference Data:
```
"gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta"
"gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai"
"gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict"
```


## Example Output
The output file is a .bam file.
The output index file is a .bam.bai file.
The validation output file is a .txt file.

Output Bam:
(INPUT ONCE THE WDL HAS RUN)
<br>
Output Bam Index:
(INPUT ONCE THE WDL HAS RUN)
<br>
Validation Output .txt file:
(INPUT ONCE THE WDL HAS RUN)
<br>


### Limitations and Important Notes

This task is relatively expensive to run, taking around 11 hours and costing $5 per sample. It is recommended that you run this task on a single sample before inputting a larger sample set. 

## Version
### 0.0.0 (10-10-2017)
- (Outdated) Initial WDL, created with WDL dev.
### 1.0.0 (07-25-2024)
- Updated to WDL version 1.0
- Utilizes the latest docker images (broadinstitute/gatk:4.6.0.0)
- Cut out unnecessary/repetitive code
- Fixed bug that caused the ValidateBam task to fail
