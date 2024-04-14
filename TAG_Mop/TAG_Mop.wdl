version 1.0

workflow TAG_Mop{
        input{
            String namespace = "broadtagteam"
            String workspaceName
            String mopDocker = "us.gcr.io/tag-team-160914/neovax-parsley:2.2.1.0"
            Boolean runMop
        }

    call rmSysfiles {
            input:
                namespace = namespace,
                workspaceName = workspaceName,
                mopDocker = mopDocker
        }

        if (runMop){
            call mop {
                input:
                    namespace = namespace,
                    workspaceName = workspaceName,
                    mopDocker = mopDocker,
                    sysfiles = rmSysfiles.deleted_sys_files

            }
        }


        output{
            Int num_deleted_sys_files = rmSysfiles.deleted_sys_files
            File deleted_sys_files = rmSysfiles.sys_files_to_delete
            Int? num_mopped_files = mop.num_of_files_to_mop
            File? mopped_files = mop.mopped_files
        }

        meta {
            author: "Yueyao Gao"
            email: "gaoyueya@broadinstitute.org"
            description: "TAG Mop contains three sub-workflows: rmSysfiles and mop. rmSysfiles removes system files that were generated from submissions from a Terra workspace. mop runs the FISS Mop function. Suggest to run after cleanupFailedSubmission.wdl"
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
            import subprocess

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
            with open('sys_files_to_delete.txt', 'w') as f:
                for file in sys_files_to_delete:
                    f.write(file + '\n')

            if len(sys_files_to_delete) == 0:
                print("No system files to delete")
            else:
                for pattern in set([i.split('/')[-1] for i in sys_files_to_delete]):
                    subprocess.run(['gsutil', '-m', 'rm', f'gs://{bucket_name}/**/{pattern}'])

            CODE

        >>>
        output{
            Int deleted_sys_files = read_int("num_of_sys_files_to_delete.txt")
            File sys_files_to_delete = "sys_files_to_delete.txt"
        }
        runtime {
            docker: mopDocker
            memory: "32 GiB"
            cpu: 8
        }
    }

    task mop {
        input{
            String namespace
            String workspaceName
            String mopDocker
            Int sysfiles
        }
        command <<<
            source activate NeoVax-Input-Parser
            # The number of system files that were deleted
            echo "System Files Deleted: ~{sysfiles}"
            # Dry run Mop
            fissfc mop -w ~{workspaceName} -p ~{namespace} --dry-run > mop_dry_run.txt
            echo Files to mop:" $(cat mop_dry_run.txt | wc -l)"
            cat mop_dry_run.txt | wc -l > num_of_files_to_mop.txt

            # Mop
            if [ $(cat mop_dry_run.txt | wc -l) -eq 0 ]; then
                echo "No files to mop"
            else
                fissfc mop -w ~{workspaceName} -p ~{namespace}
            fi
        >>>
        output{
            Int num_of_files_to_mop = read_int("num_of_files_to_mop.txt")
            File mopped_files = "mop_dry_run.txt"
        }
        runtime {
            docker: mopDocker
            memory: "32 GiB"
            cpu: 8
        }
    }
