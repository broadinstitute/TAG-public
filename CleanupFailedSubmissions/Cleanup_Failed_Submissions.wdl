version 1.0

workflow Cleanup_Failed_Submissions {
    # This script is based on CleanupIntermediate.wdl from long-read-pipeline team
    # I have added a step to get the workspace bucket name and the failed submissions.
    # Ironicaly, this generates intermeidate files too, but they are tiny.
    meta {
        description: "A workflow to clean up storage that associated with failed submission. Use at your own risk."
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
    }
    parameter_meta {
        namespace: "project to which workspace belongs (str)"
        workspace: "Terra workspace name (str)"
    }

    input {
        String namespace
        String workspace
        Boolean remove_partially_fail = false
    }

    call GetWorkspaceInfo {
        input:
            namespace = namespace,
            workspace = workspace,
            remove_partially_fail = remove_partially_fail

    }

    scatter (sid in GetWorkspaceInfo.failed_submissions) {
        call CleanupAFolder {
            input:
                bucket_name = GetWorkspaceInfo.workspace_bucket,
                submission_id = sid
        }
    }

output {
        Int num_failed_submissions = length(GetWorkspaceInfo.failed_submissions)
    }
}

task GetWorkspaceInfo {
    input {
        String namespace
        String workspace
        Boolean remove_partially_fail
    }
    command <<<
        source activate NeoVax-Input-Parser
        python3 <<CODE

        import firecloud.api as fapi

        namespace = "~{namespace}"
        workspace = "~{workspace}"

        if '~{remove_partially_fail}' == 'true':
            with open('failed_submissions.txt','w') as file:
                for submission in fapi.list_submissions(namespace, workspace).json():
                    if 'Failed' in submission['workflowStatuses'].keys() or 'Aborted' in submission['workflowStatuses'].keys():
                        file.write(submission['submissionId'] + '\n')
        else:
            with open('failed_submissions.txt','w') as file:
                    for submission in fapi.list_submissions(namespace, workspace).json():
                        if 'Failed' in submission['workflowStatuses'].keys() or 'Aborted' in submission['workflowStatuses'].keys():
                            if 'Failed' and 'Succeeded' not in submission['workflowStatuses'].keys():
                                file.write(submission['submissionId'] + '\n')

        with open("workspace_bucket.txt", "w") as file:
            file.write(fapi.get_workspace(namespace, workspace).json()['workspace']['bucketName'])

        CODE
        >>>
    runtime {
       docker: "us.gcr.io/tag-team-160914/neovax-parsley:2.2.1.0"
       preemptible: 0
    }
    output {
        Array[String] failed_submissions = read_lines("failed_submissions.txt")
        String workspace_bucket = read_string("workspace_bucket.txt")

    }
}

task CleanupAFolder {
    input {
        String bucket_name
        String submission_id
    }

    command <<<
        # Older version of Terra does not have the submission folder
        if gsutil -q ls gs://~{bucket_name}/submissions/~{submission_id} >/dev/null 2>&1; then
            timeout 23h gsutil -q rm -rf gs://~{bucket_name}/submissions/~{submission_id} || echo "Timed out. Please try again."
        elif gsutil -q ls gs://~{bucket_name}/~{submission_id} >/dev/null 2>&1; then
            timeout 23h gsutil -q rm -rf gs://~{bucket_name}/~{submission_id} || echo "Timed out. Please try again."
        else
            echo "File path does not exist."
        fi
    >>>

    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 10 HDD"
        preemptible_tries:     1
        max_retries:           1
        docker:"us.gcr.io/google.com/cloudsdktool/google-cloud-cli:alpine"
    }
}
