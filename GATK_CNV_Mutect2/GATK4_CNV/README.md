# GATK_CNV
Workflow for calling somatic copy number variants as tumor-normal or tumor-only mode. Supports both WGS and WES.

This WDL is originally from [CNV_Somatic_Panel_Workflow](https://github.com/gatk-workflows/gatk4-somatic-cnvs) and 
is modified to output the number of segments, amplifications, and deletions in a Terra workspace table.


## Notes:
* The intervals argument is required for both WGS and WES workflows and accepts formats compatible with the
   GATK -L argument (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
   These intervals will be padded on both sides by the amount specified by padding (default 250)
   and split into bins of length specified by bin_length (default 1000; specify 0 to skip binning,
   e.g., for WES).  For WGS, the intervals should simply cover the autosomal chromosomes (sex chromosomes may be
   included, but care should be taken to 1) avoid creating panels of mixed sex, and 2) denoise case samples only
   with panels containing only individuals of the same sex as the case samples).

* Intervals can be blacklisted from coverage collection and all downstream steps by using the blacklist_intervals
   argument, which accepts formats compatible with the GATK -XL argument
   (see https://gatkforums.broadinstitute.org/gatk/discussion/11009/intervals-and-interval-lists).
   This may be useful for excluding centromeric regions, etc. from analysis.  Alternatively, these regions may
   be manually filtered from the final callset.

* A reasonable blacklist for excluded intervals (-XL) can be found at:
   hg19: gs://gatk-best-practices/somatic-b37/CNV_and_centromere_blacklist.hg19.list
   hg38: gs://gatk-best-practices/somatic-hg38/CNV_and_centromere_blacklist.hg38liftover.list (untested)

* The sites file (common_sites) should be a Picard or GATK-style interval list.  This is a list of sites
   of known variation at which allelic counts will be collected for use in modeling minor-allele fractions.

## Example invocation
```angular2html
java -jar cromwell.jar run GATK_CNV.wdl -i my_parameters.json
```

