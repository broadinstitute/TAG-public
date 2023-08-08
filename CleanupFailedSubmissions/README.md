# Cleanup Failed Submissions Workflow

## Description

This workflow is based on the `CleanupIntermediate.wdl` script from the long-read-pipeline team, with an additional step to retrieve the workspace bucket name and failed submissions. This workflow is designed to clean up storage associated with failed submissions. Please use it cautiously and at your own risk.

## Workflow Usage

This workflow requires the input of the `namespace` and `workspace` where the failed submissions need to be cleaned up.


- `namespace`: Project to which the workspace belongs 
- `workspace`: Terra workspace name

### Tasks

1. **GetWorkspaceInfo:** This task retrieves the workspace bucket name and identifies failed submissions using the FireCloud API. It generates the `failed_submissions.txt` file containing submission IDs and the `workspace_bucket.txt` file containing the workspace bucket name.

2. **CleanupAFolder:** This task cleans up a specified folder associated with a submission. It removes the folder from the corresponding bucket using `gsutil`.

### Workflow Steps

1. Call the **GetWorkspaceInfo** task to retrieve failed submissions and workspace bucket name.
2. Scatter over each failed submission ID.
3. For each failed submission, call the **CleanupAFolder** task to remove the submission's associated folder from the workspace bucket.

## Instructions

1. Copy this workflow into your Terra workspace.
2. Provide the `namespace` and `workspace` as inputs to the workflow.
3. Execute the workflow, and it will clean up storage associated with failed submissions.

Please use this workflow carefully, as it involves the deletion of data.

## Dependencies
This workflow is intended to be used by TAG team member internally. If you can't access the following docker image,
please create a docker image with Firecloud API installed and replace the docker image in the workflow.
- Docker image: `us.gcr.io/tag-team-160914/neovax-parsley:2.2.1.0`

---
For questions or concerns, contact the author at tag@broadinstitute.org.
