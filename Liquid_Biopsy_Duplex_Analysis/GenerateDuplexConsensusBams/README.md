# Generate Duplex Consensus Bams
During the library construction of liquid biopsy samples, Unique Molecular Identifiers (UMIs) are added to both ends of reads for sequencing to obtain deeper coverage and better quality reads. This workflow performs UMI-aware de-duplication and generates analysis-ready consensus BAMs by taking advantage of the UMI information.

## Note:
This workflow was copied from [GP-TAG/GenerateDuplexConsensusBams](https://portal.firecloud.org/?return=terra#methods/GP-TAG/GenerateDuplexConsensusBams/19)
snapshot 19