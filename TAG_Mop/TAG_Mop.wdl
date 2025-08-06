version 1.0

workflow CleanupWithOptionalMop {
    input {
        String namespace
        String workspaceName
        String mopDocker = "us.gcr.io/tag-team-160914/neovax-parsley:2.2.1.0"
        String? allowed_submitters
        Boolean delete_sys_files
    }

    # Optional mop calls
    if (defined(allowed_submitters)) {
        call FilterSubmissionIdsBySubmitter {
            input:
                namespace = namespace,
                workspaceName = workspaceName,
                allowed_submitters = allowed_submitters,
                mopDocker = mopDocker
        }

        call mop as mop_with_submitter {
            input:
                namespace = namespace,
                workspaceName = workspaceName,
                mopDocker = mopDocker,
                submission_ids_file = FilterSubmissionIdsBySubmitter.filtered_submission_ids
        }
    }

    if (!defined(allowed_submitters)) {
        call mop as mop_without_submitter {
            input:
                namespace = namespace,
                workspaceName = workspaceName,
                mopDocker = mopDocker
        }
    }

    # Always define a dummy fallback trigger
    String dummy_trigger = "noop"

    if (delete_sys_files && defined(allowed_submitters)) {
        call rmSysfiles as rm_after_mop_with_submitter {
            input:
                namespace = namespace,
                workspaceName = workspaceName,
                mopDocker = mopDocker,
                trigger = mop_with_submitter.total_size_to_mop
        }
    }

    if (delete_sys_files && !defined(allowed_submitters)) {
        call rmSysfiles as rm_after_mop_without_submitter {
            input:
                namespace = namespace,
                workspaceName = workspaceName,
                mopDocker = mopDocker,
                trigger = mop_without_submitter.total_size_to_mop
        }
    }

    output {
        Int deleted_sysfiles = if delete_sys_files then select_first(select_all([
            rm_after_mop_with_submitter.deleted_sys_files,
            rm_after_mop_without_submitter.deleted_sys_files
        ])) else 0

        File? mopped_files = select_first([
            mop_with_submitter.mopped_files,
            mop_without_submitter.mopped_files
        ])

        String? total_size_to_mop = select_first([
            mop_with_submitter.total_size_to_mop,
            mop_without_submitter.total_size_to_mop
        ])
    }


    meta {
        author: "Yueyao Gao"
        email: "gaoyueya@broadinstitute.org"
        description: "TAG Mop optionally runs mop and/or removes sysfiles after mop."
    }

 }




    task rmSysfiles {
        input{
            String namespace
            String workspaceName
            String mopDocker
            Int memory = 16
            Int cpu = 2
            String? trigger # dummy input to inforce order
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
            storage_client = storage.Client(namespace)
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
            memory: memory + " GiB"
            cpu: cpu
        }
    }


task FilterSubmissionIdsBySubmitter {
    input {
        String namespace
        String workspaceName
        String? allowed_submitters
        String mopDocker
    }

    command <<<
        source activate NeoVax-Input-Parser
        export NAMESPACE="~{namespace}"
        export WORKSPACE="~{workspaceName}"
        export ALLOWED_SUBMITTERS="~{select_first([allowed_submitters, ""])}"

        python3 <<CODE
import os
import firecloud.api as fapi

namespace = os.environ["NAMESPACE"]
workspace = os.environ["WORKSPACE"]
allowed_raw = os.environ.get("ALLOWED_SUBMITTERS", "")
allowed_submitters = [s.strip() for s in allowed_raw.split(",") if s.strip()]

submissions = fapi.list_submissions(namespace, workspace).json()

filtered_ids = [
    s["submissionId"]
    for s in submissions
    if not allowed_submitters or s.get("submitter", "") in allowed_submitters
]

with open("filtered_submission_ids.txt", "w") as f:
    for sid in filtered_ids:
        f.write(sid + "\\n")
CODE
    >>>

    output {
        File filtered_submission_ids = "filtered_submission_ids.txt"
    }

    runtime {
        docker: mopDocker
        memory: "2 GiB"
        cpu: 1
    }
}



    task mop {
        input {
            String namespace
            String workspaceName
            String mopDocker
            File? submission_ids_file
        }

        command <<<
            source activate NeoVax-Input-Parser

            if [[ -f "~{submission_ids_file}" ]]; then
                SUB_IDS=$(paste -sd " " ~{submission_ids_file})
                fissfc mop -w ~{workspaceName} -p ~{namespace} --dry-run --submission-ids $SUB_IDS > mop_dry_run.txt
            else
                fissfc mop -w ~{workspaceName} -p ~{namespace} --dry-run > mop_dry_run.txt
            fi

            echo Files to mop: $(cat mop_dry_run.txt | wc -l)
            cat mop_dry_run.txt | grep "gs:" | wc -l > num_of_files_to_mop.txt
            cat mop_dry_run.txt | grep Total | awk '{print $3 $4}' > total_size_to_mop.txt

            if [ $(cat mop_dry_run.txt | wc -l) -eq 0 ]; then
                echo "No files to mop"
            else
                if [[ -f "~{submission_ids_file}" ]]; then
                    fissfc --yes mop -w ~{workspaceName} -p ~{namespace} --submission-ids $SUB_IDS
                else
                    fissfc --yes mop -w ~{workspaceName} -p ~{namespace}
                fi
            fi
        >>>

        output {
            File mopped_files = "mop_dry_run.txt"
            Int num_of_files_to_mop = read_int("num_of_files_to_mop.txt")
            String total_size_to_mop = read_string("total_size_to_mop.txt")
        }

        runtime {
            docker: mopDocker
            memory: "4 GiB"
            cpu: 1
        }
    }