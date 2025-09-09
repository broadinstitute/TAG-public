# BrightSeq Somatic WES CNV Analysis
For BrightSeq WES, TAG is managing a workspace that will automatically launch the GATK4 CNV workflow on tumor-normal pairs as they are delivered. The CNVSomaticPairWorkflow_BCH_Wrapper version of the CNV workflow also runs QuicViz and evaluates the CNV for clinical PASS/FAIL criteria. 

The DeliverCNVOutputs WDL is set to copy the clinical deliverables into the collaborator's workspace and will be run by the PM after they have delivered the short variant data to the collaborator.