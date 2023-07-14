# SvABA Somatic Structural Variant Caller
SvABA detectes somatic rearrangements and indels by first assembling contigs from reads that diverge from the reference sequence, then aligning those extended contigs back to the reference to detect SV events.

This workflow calls somatic SVs with SvABA v0.2.1 for exome, custom panel, or whole-genome BAMs.  SvABA can perform either tumor-only or tumor-normal analysis.  For the tumor-only mode, the tool requires an unmatched normal.  The workflow has a list of known germline variants and artifacts for reducing false positive SV calls in tumor samples.

### Citation
Wala JA et al. SvABA: genome-wide detection of structural variants and indels by local assembly. Genome Res. 2018 Apr;28(4):581-591. doi: 10.1101/gr.221028.117. Epub 2018 Mar 13. PMID: 29535149; PMCID: PMC5880247.