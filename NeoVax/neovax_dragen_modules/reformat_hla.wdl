version 1.0

task reformat_hla {
	input {
		String patient_id
		File hla_file
	}
	command <<<
python <<CODE
allele_type = ["A", "B", "C"]
allele_start = "HLA-"
output_file = "~{patient_id}.reformated_hla_alleles.txt"
hla_alleles = open("~{hla_file}", 'r').readlines()[0].strip().split('\t')[-1].split(' ')
with open(output_file, 'w') as writer:
	for allele in hla_alleles:
		if allele[0] in allele_type:
			formatted_allele =  allele_start + allele[0] + allele[2:]
			writer.write(formatted_allele + '\n')
CODE
	>>>

	output {
		File reformatted_hla = "~{patient_id}.reformated_hla_alleles.txt"
	}

	runtime {
		docker: "python:3"
		memory: "2 GB"
		disks: "local-disk 10 HDD"
		preemptible: 1
		maxRetries: 1
	}
}

workflow reformat_hla_workflow {
	call reformat_hla
}