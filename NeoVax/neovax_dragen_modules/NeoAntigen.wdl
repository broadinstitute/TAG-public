version 1.0

workflow Neoantigen_Pipeline {
	input {
		String id
		String caseName
		File hg19DBTarBall
		File mafFile
		File geneExpr
		File duplicate_genes
		File gene_id_name_mapping
	}

	call neoantigen {
	  input:
		id = id,
		caseName = caseName,
		hg19DBTarBall = hg19DBTarBall,
		mafFile = mafFile,
		geneExpr = geneExpr
	}
	call generate_peptides {
	  input:
		pairName = id,
		neoorfs = neoantigen.neoorfFile,
		snps = neoantigen.outSNP_for_ms_predictions,
		indels = neoantigen.outINDEL_for_ms_predictions,
		expression = geneExpr,
		duplicate_genes = duplicate_genes,
		gene_id_name_mapping = gene_id_name_mapping
	}
	output {
		File outINDEL_for_ms_predictions = neoantigen.outINDEL_for_ms_predictions
		File outSNP_for_ms_predictions = neoantigen.outSNP_for_ms_predictions
		File combinedFullpeptideFile = neoantigen.combinedFullpeptideFile
		File neoorfFile = neoantigen.neoorfFile
		File snp_indel_neoantigens = generate_peptides.snp_indel_neoantigens
		File neoorf_neoantigens = generate_peptides.neoorf_neoantigens
		File nonsense_vars_file = generate_peptides.nonsense_vars_file
	}
}


task generate_peptides {
	input {
		File snps
		File indels
		File neoorfs
		File? expression
		String pairName
		File duplicate_genes
		File gene_id_name_mapping

		Int preemptible = 2
		Int maxRetries =1
		Int memoryGB = 3
		Int diskGB = 50
		String? docker_override
	}
	String nonsense_var_table = pairName + ".nonsense_vars.txt"

	Boolean expression_provided = defined(expression)

	command <<<
		set -x

		if [ ~{expression_provided} = true ]; then

		python<<CODE
		import pandas as pd
		rsem = pd.read_csv("~{expression}",sep='\t')
		gene_id_name_mapping = pd.read_csv("~{gene_id_name_mapping}",sep='\t',names=['gene_id','Hugo_Symbol'])
		merged = rsem.merge(gene_id_name_mapping,on='gene_id',how='inner')
		merged.to_csv("~{pairName}_tpm_hugo_symbol.txt",sep='\t',index=False)
		CODE

		grep -f ~{duplicate_genes} "~{pairName}_tpm_hugo_symbol.txt"

		fi

		touch ~{nonsense_var_table}

		Rscript /format.R \
		--snps ~{snps} \
		--indels ~{indels} \
		--identifier ~{pairName} \
		--gen_path /generate_peptides.R \
		--templates /template_cols \
		--nonsense_file ~{nonsense_var_table} \
		--output_dir .

		for length in 8 9 10 11; do
			python /parse_updn.py ~{pairName}.peps.txt $length
		done

		python /parse_neoorf_updn.py ~{neoorfs} ~{pairName}_neoorfs

		head -n 1 ~{pairName}.peps9_parsed.tsv > snp_indel_header.txt
		head -n 1 ~{pairName}_neoorfs9.txt > neoorf_header.txt
		cat ~{pairName}.peps*_parsed.tsv | grep -vf snp_indel_header.txt > snp_indel_noheader.txt
		cat ~{pairName}_neoorfs*.txt | grep -vf neoorf_header.txt > neoorfs_noheader.txt
		cat snp_indel_header.txt snp_indel_noheader.txt > snp_indel.txt
		cat neoorf_header.txt neoorfs_noheader.txt > neoorfs.txt

		if [ ~{expression_provided} = true ]; then

		python<<CODE
		import pandas as pd
		snp_indel = pd.read_csv('snp_indel.txt',sep='\t')
		neoorf = pd.read_csv('neoorfs.txt',sep='\t')
		expr = pd.read_csv("~{pairName}_tpm_hugo_symbol.txt",sep='\t')
		snp_indel_expr = snp_indel.merge(expr,on='Hugo_Symbol',how='inner')[['Hugo_Symbol','pep_wt','pep','ctex_up','ctex_dn','gene_id','TPM','sample','Entrez_Gene_Id','Chromosome','Start_Position','End_Position','Strand','Variant_Classification','Variant_Type','Annotation_Transcript','Genome_Change','cDNA_Change','Codon_Change','Protein_Change']]
		neoorf_expr = neoorf.merge(expr,on='Hugo_Symbol',how='inner')[['Hugo_Symbol','pep_wt','pep','ctex_up','ctex_dn','gene_id','TPM','sample','Entrez_Gene_Id','Chromosome','Start_Position','End_Position','Strand','Variant_Classification','Variant_Type','Annotation_Transcript','Genome_Change','cDNA_Change','Codon_Change','Protein_Change']]
		snp_indel_expr.to_csv('~{pairName}.snp_indel.tsv',sep='\t',index=False)
		neoorf_expr.to_csv('~{pairName}.neoorfs.tsv',sep='\t',index=False)
		CODE

		else

			mv snp_indel.txt "~{pairName}.snp_indel.tsv"
			mv neoorfs.txt "~{pairName}.neoorfs.tsv"

		fi


	>>>

	output {
		Array[File] all_outs = glob("*")
		File snp_indel_neoantigens = "~{pairName}.snp_indel.tsv"
		File neoorf_neoantigens = "~{pairName}.neoorfs.tsv"
        File nonsense_vars_file = "~{nonsense_var_table}"
	}


	runtime {
		docker: select_first([docker_override, "us.gcr.io/tag-team-160914/neovax-tag-format-ms_intrinsic:v2"])
		memory: memoryGB+" GB"
		disks: "local-disk " + diskGB + " HDD"
		preemptible: preemptible
		maxRetries: maxRetries
	}
}

task neoantigen {
	input {
		# RUNTIME INPUT PARAMS
		Int preemptible = 1
		Int memoryGB = 60
		Int diskGB_buffer = 20

		# TASK INPUT PARAMS
		File hg19DBTarBall
		File mafFile
		File geneExpr
		String id
		String caseName
		String? docker_override
	}
	# COMPUTE DISK SIZE
	Int diskGB = ceil(size(hg19DBTarBall, 'G') * 4 + size(mafFile, 'G') + diskGB_buffer)


	command <<<
		set -euxo pipefail
		# correct [Start|End]_Position to [Start|End]_position
		grep ^Hugo_Symbol ~{mafFile} | tr "\t" "\n" > ~{id}.column_names.txt
		python /root/neoantigen/src/subset_maf_columns.py --correct-maf-dialect \
		~{id}.column_names.txt ~{mafFile} ~{id}.corrected.maf

		mutFile="~{id}.snp.maf"
		indelFile="~{id}.indel.maf"
		cat <(grep ^Hugo ~{id}.corrected.maf) \
			<(awk -F "\t" '{if($10 != "INS" && $10 != "DEL"){print $0}}' OFS="\t" ~{id}.corrected.maf) > $mutFile
		cat <(grep ^Hugo ~{id}.corrected.maf) \
			<(awk -F "\t" '{if($10 == "INS" || $10 == "DEL"){print $0}}' OFS="\t" ~{id}.corrected.maf) > $indelFile

		HG19_DB_DIR_NAME=$(tar tvf ~{hg19DBTarBall} | egrep -o "[^ ]+/$")
		tar -xvf ~{hg19DBTarBall}

		# modified get_mut_protein_gencode_v1.pl
		/root/neoantigen/src/get_mut_protein_gencode_v1.pl $mutFile infile transcript_gencode_v19 ~{id} `pwd`/$HG19_DB_DIR_NAME > ~{id}.mutFile.out
		/root/neoantigen/src/get_mut_protein_gencode_v1.pl $indelFile infile transcript_gencode_v19 ~{id}.indel `pwd`/$HG19_DB_DIR_NAME > ~{id}.indelFile.out

		cat ~{id}.fullpeptide.fa ~{id}.indel.fullpeptide.fa > ~{id}.combined.fullpeptide.fa

		echo "#### get_neoorf_table ####"
		/root/neoantigen/src/get_neoorf_table.pl ~{id}.combined.fullpeptide.fa /root/neoantigen/wt.peptide.9mers.txt /root/neoantigen/wt.peptide.10mers.txt

		/root/neoantigen/src/process_fullpeptide_for_ms_predictions.pl $mutFile ~{id}.pep81.fa ~{caseName} > ~{id}.hg19.snp.pep81.for_ms_predictions.txt
		/root/neoantigen/src/process_fullpeptide_for_ms_predictions.pl $indelFile ~{id}.indel.pep81.fa ~{caseName} > ~{id}.hg19.indel.pep81.for_ms_predictions.txt

	>>>

	runtime {
		docker: select_first([docker_override, "us.gcr.io/tag-team-160914/neovax-tag-neoantigen:v3"])
		memory: "~{memoryGB} GB"
		disks: "local-disk ~{diskGB} HDD"
		preemptible: "~{preemptible}"
	}

	output {
		File combinedFullpeptideFile = "~{id}.combined.fullpeptide.fa"
		File neoorfFile = "~{id}.combined.neoorf.table.txt"
		File outSNP_for_ms_predictions = "~{id}.hg19.snp.pep81.for_ms_predictions.txt"
		File outINDEL_for_ms_predictions = "~{id}.hg19.indel.pep81.for_ms_predictions.txt"
	}
}
