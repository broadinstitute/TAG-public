version 1.0

workflow TAG_Mop {
    input {
        String namespace = "broadtagteam"
        String workspaceName
        String mopDocker = "us.gcr.io/tag-team-160914/neovax-parsley:2.2.1.0"
        Boolean removeFailedSubmissions
        Boolean runMop
        Boolean remove_partially_fail = false
    }

    scatter (task_index in range(3)) {
        if (task_index == 0) {
            call rmSysfiles {
                input:
                    namespace = namespace,
                    workspaceName = workspaceName,
                    mopDocker = mopDocker
            }
        }
        if (removeFailedSubmissions && task_index == 1) {
            call GetFailedSubmissions {
                input:
                    namespace = namespace,
                    workspaceName = workspaceName,
                    mopDocker = mopDocker,
                    remove_partially_fail = removeFailedSubmissions
            }
            scatter (sid in GetFailedSubmissions.failed_submissions) {
                call CleanupAFolder {
                    input:
                        bucket_name = GetFailedSubmissions.workspace_bucket,
                        submission_id = sid
                }
            }
        }
        if (runMop && task_index == 2) {
            call mop {
                input:
                    workspaceName = workspaceName,
                    mopDocker = mopDocker
            }
        }
    }

    output {
        Int deleted_sys_files = rmSysfiles.deleted_sys_files
    }

    meta {
        author: "Yueyao Gao"
        email: "gaoyueya@broadinstitute.org"
        description: "TAG Mop contains three sub-workflows: rmSysfiles, removeFailedSubmission, and mop. rmSysfiles removes system files that were generated from submissions from a Terra workspace. mop runs the Mop pipeline."
    }
}

    task rmSysfiles {
        input{
            String namespace
            String workspaceName
            String mopDocker
        }
        command <<<
            source activate NeoVax-Input-Parser
            python <<CODE
            from google.cloud import storage
            import firecloud.api as fapi

            namespace = "~{namespace}"
            workspaceName = "~{workspaceName}"
            bucket_name = fapi.get_workspace(namespace, workspaceName).json()['workspace']['bucketName']

            # Collect the system files to delete
            storage_client = storage.Client()
            blobs = storage_client.list_blobs(bucket_name, projection='full')
            patterns_to_remove = ["stdout.log", "stderr.log", "localization.sh", "gcs_transfer.sh", "/stdout","/stderr","/rc","-rc.txt",'/memory_retry_rc','/output','/script','/exec.sh']
            sys_files_to_delete = []
            for blob in blobs:
                for pattern in patterns_to_remove:
                    if blob.name.endswith(pattern):
                        sys_files_to_delete.append(f"gs://{bucket_name}/{blob.name}")

            # Output the number of system files to delete
            with open('num_of_sys_files_to_delete.txt', 'w') as f:
                f.write(str(len(sys_files_to_delete)))
            print(f"System Files to Delete in {namespace}/{workspaceName}: ", len(sys_files_to_delete))

            for pattern in set([i.split('/')[-1] for i in sys_files_to_delete], desc="Deleting System Files", unit="pattern"):
                subprocess.run(['gsutil', '-m', 'rm', f'gs://{bucket_name}/**/{pattern}'])

            CODE

        >>>
        output{
            Int deleted_sys_files = read_int("num_of_sys_files_to_delete.txt")
        }
        runtime {
            docker: mopDocker
            memory: "32 GiB"
            cpu: 8
        }
    }

    task GetFailedSubmissions {
        input {
            String namespace
            String workspaceName
            Boolean remove_partially_fail
            String mopDocker
        }
        command <<<
            source activate NeoVax-Input-Parser
            python3 <<CODE

            import firecloud.api as fapi

            namespace = "~{namespace}"
            workspace = "~{workspaceName}"

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
           docker: mopDocker
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
        timeout 23h gsutil -q rm -rf gs://~{bucket_name}/submissions/~{submission_id} || echo "Timed out. Please try again."
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


    task mop {
        input{
            String namespace
            String workspaceName
            String mopDocker
        }
        command <<<
            source activate NeoVax-Input-Parser
            python <<CODE
            import firecloud.api as fapi
            import subprocess

            namespace = "~{namespace}"
            workspaceName = "~{workspaceName}"

            # Run fissfc Mop to remove data that not presented in the data model
            subprocess.run(['fissfc', 'mop', '-w', workspaceName, '-p', namespace])

            CODE

        >>>
        runtime {
            docker: mopDocker
            memory: "32 GiB"
            cpu: 8
        }
    }
