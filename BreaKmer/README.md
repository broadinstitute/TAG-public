# BreaKmer Structural Variant Caller
BreaKmer calls structural variants with annotated breakpoints by assembling contigs from reads which contain unique k-mer sequences not found in the reference and then re-aligning those extended contigs to the reference. Because BreaKmer will attempt to assemble contigs from all heavily soft-clipped reads, the workflow first downsamples and re-aligns the input BAMs while clipping overhangs to reduce the runtime of the tool. 

### Citation
Abo RP, Ducar M, Garcia EP, Thorner AR, Rojas-Rudilla V, Lin L, Sholl LM, Hahn WC, Meyerson M, Lindeman NI, Van Hummelen P, MacConaill LE. BreaKmer: detection of structural variation in targeted massively parallel sequencing data using kmers. Nucleic Acids Res. 2015 Feb 18;43(3):e19. doi: 10.1093/nar/gku1211. Epub 2014 Nov 26. PMID: 25428359; PMCID: PMC4330340.