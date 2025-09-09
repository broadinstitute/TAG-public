# Signature Profiling WDL for CODEC Mutlist Output
CODEC pipeline: SingleSampleCODEC pipeline provides text files for discovered mutations. 

This workflow summarizes and plots mutation spectrums in 96 trinucleotide contexts and generate Mutation Matrix that will be later used to subtract SBS(Single Base Substitution) signatures from SNVs with database reference from COSMIC(https://cancer.sanger.ac.uk/signatures/).

The Signature Profiling tool is from https://github.com/AlexandrovLab/SigProfilerAssignment and has been implanted to the docker image.

The output of this WDL includes: 
1) SpectrumPlots
2) MutationMetrics
3) SignatureCount
4) SignatureProportionPDF
5) SignatureStackedPlot
6) TMBPlot
7) DecomposedSignatureProbabilities


### Citation
Díaz-Gay et al. 2023 Bioinformatics and Tate et al. 2019 Nucleic Acids Research