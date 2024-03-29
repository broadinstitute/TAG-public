# Mutect2-GATK4
This WDL workflow runs GATK4 Mutect2 on a single tumor-normal pair or on a single tumor sample, and performs additional filtering and functional annotation tasks.

## Main requirements/expectations :

* One analysis-ready BAM file (and its index) for each sample

## Description of inputs:

### Runtime

* gatk_docker, oncotator_docker: docker images to use for GATK 4 Mutect2 and for Oncotator
* tag_docker: docker images for TAG's add-on analyses
* preemptible_attempts: how many preemptions to tolerate before switching to a non-preemptible machine (on Google)
* *_override: (optional) local file or Google bucket path to scripts to be used instead of the scripts in the docker image

### Workflow options

* intervals: genomic intervals (will be used for scatter)
* scatter_count: number of parallel jobs to generate when scattering over intervals
* artifact_modes: types of artifacts to consider in the orientation bias filter (optional)
* m2_extra_args, m2_extra_filtering_args: additional arguments for Mutect2 calling and filtering (optional)
* split_intervals_extra_args: additional arguments for splitting intervals before scattering (optional)
* run_orientation_bias_filter: if true, run the orientation bias filter post-processing step (optional, false by default)
* run_oncotator: if true, annotate the M2 VCFs using oncotator (optional, false by default)

### Primary inputs

* ref_fasta, ref_fai, ref_dict: reference genome, index, and dictionary
* tumor_bam, tumor_bam_index: BAM and index for the tumor sample
* normal_bam, normal_bam_index: BAM and index for the normal sample

### Primary resources (optional but strongly recommended)

* pon, pon_index: optional panel of normals in VCF format containing probable technical artifacts (false positves)
* gnomad, gnomad_index: optional database of known germline variants (see http://gnomad.broadinstitute.org/downloads)
* variants_for_contamination, variants_for_contamination_index: VCF of common variants with allele frequencies for calculating contamination* 

### Secondary resources (for optional tasks)

* onco_ds_tar_gz, default_config_file: Oncotator datasources and config file
* sequencing_center, sequence_source: metadata for Oncotator
* filter_oncotator_maf: Whether the MAF generated by oncotator should have the filtered variants removed (true by default)

### TAG's modification

* HaplotypeCaller: task runs HaplotypeCaller on normal bam
* CallableLoci: task uses GATK3 CallableLoci to compute the number of somatically callable bases and 3-base mutational contexts
* MutationalBurden: task reads MAF and computes both coding and non-coding mutational burdens (# of mutations / callable bases)
* LegoPlot: task generates lego plots

