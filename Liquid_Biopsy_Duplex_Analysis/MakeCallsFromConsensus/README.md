# MakeCallsFromConsensus
This workflow calls SNVs (INDEL calling is unsupported in this version) 
from duplex consensus BAMs. If you haven't created the consensus BAMs, 
please do so by checking GenerateDuplexConsensusBams workflow.
You can run the workflow in either tumor-normal 
(if normal samples are available) or tumor-only mode.

## Note:
This workflow wdl was copied from [GP-TAG/MakeCallsFromConsensus](
https://portal.firecloud.org/?return=terra#methods/GP-TAG/MakeCallsFromConsensus/47/wdl
) snapshot 47