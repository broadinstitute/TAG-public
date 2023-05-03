# Griffin WDL
This directory provides a set of WDL scripts to streamline the execution of [Griffin](https://github.com/adoebley/Griffin) 
on ULP-WGS data for Metastatic Breast Cancer (MBC) subtype classification. These WDL scripts will allow you to run the 
Griffin workflow on Terra. This pipeline is composed of the following steps:
* **GC correction**: 
Calculate the GC bias for a given set of bam files
* **Nucleosome profiling**:
Run nucleosome profiling for a given set of site lists and a given set of bam files
* **ER classifier**: 
Run ER classifier on features generated from the nucleosome profiling steps

You can find an example Terra workspace that uses Griffin WDL in the [broadtagteam/tag1592_Griffin_MBC_ULP_WGS](https://app.terra.bio/#workspaces/broadtagteam/tag1592_Griffin_MBC_ULP_WGS).

## IMPORTANT NOTE
* Griffin only support hg38 reference genome. If your data is in hg19, you will need to realign your bam files to hg38.
* The provided input json is a ready to use example json template of the Griffin workflow. It is meant to run on a demo 
bam to test the functionality of griffin_GC_correction and griffin_nucleosome_profiling.
* You will need to modify the input json to run on your own data. Please check [broadtagteam/tag1592_Griffin_MBC_ULP_WGS](https://app.terra.bio/#workspaces/broadtagteam/tag1592_Griffin_MBC_ULP_WGS) for more details about
MBC subtype classification json configuration.

