version 1.0

struct Runtime {
  Int max_retries
  Int preemptible
  Int mem
  Int cpu
  Int boot_disk_size
  Int initial_disk_size
  String docker
}

workflow HLAthena {
  input {
    String analysis_name

    File alleles_file
    File peptide_list
    File? second_peptide_list

    String peptide_col_name = "pep"
    Boolean aggregate_pep = false
    
    String expr_col_name = "TPM"
    Boolean exists_expr = true
    Boolean logtransform_expr = true
    Boolean exists_ctex = true

    # run NetMHCPan if executable tarball is supplied
    File? netmhc_tar_gz

    String models_path = "gs://msmodels/"
    String features_all = "features_AAPos_AAPCA_LogTPM_CNN_Kidera_Gene"
    String feature_sets = "features_AAPos_AAPCA_Kidera"

    Float assign_threshold = 2
    String assign_by_ranks_or_scores = "ranks"
    String assign_colors

    File? wild_type_peptides_tar_gz

    # workflow level runtime parameters
    Int preemptible = 2
    Int max_retries = 1
    Int additional_disk = 20
    Int boot_disk_size = 15
    Int mem = 16
    Int cpu = 1
  }

  Runtime standard_runtime = { "preemptible": preemptible,
                                "max_retries": max_retries,
                                "mem": mem,
                                "cpu": cpu,
                                "docker": "us.gcr.io/tag-public/neovax-tag-hlathena:v1",
                                "boot_disk_size": boot_disk_size,
                                "initial_disk_size": additional_disk }

  Boolean run_netmhc = if defined(netmhc_tar_gz) then true else false

  ### Parse HLA alleles
  call ParseAlleles {
    input:
      analysis_name = analysis_name,
      alleles_file = alleles_file,
      runtime_params = standard_runtime
  }
  
  ### Merge peptide lists if second list is supplied
  if (defined(second_peptide_list)) {
    call MergePeptideLists {
      input:
        analysis_name = analysis_name,
        peptide_list = peptide_list,
        second_peptide_list = select_first([second_peptide_list]),
        runtime_params = standard_runtime
    }
  }

  File input_peptides = select_first([MergePeptideLists.merged_peptide_list, peptide_list])

  ### If up/dn context is in the input file, compute and append cleavability score
  if (exists_ctex) {
    call PredictCleavability {
      input:
        analysis_name = analysis_name,
        peptide_list = input_peptides,
        peptide_col_name = peptide_col_name,
        runtime_params = standard_runtime
    }
  }

  ### Optionally aggregate peptides based on a given column
  ### (e.g. agregate all possible transcript precursors of a peptide)
  if (aggregate_pep) {
    call AggregatePeptides {
      input:
        analysis_name = analysis_name,
        peptide_list = select_first([PredictCleavability.peptides_clev, input_peptides]),
        peptide_col_name = peptide_col_name,
        exists_expr = exists_expr,
        expr_col_name = expr_col_name,
        runtime_params = standard_runtime
    }
  }

  ### If expression is in the input file, log-transform as needed
  if (exists_expr) {
    call LogTransformExpr {
      input:
        analysis_name = analysis_name,
        peptide_list = select_first([AggregatePeptides.peptides_agg, PredictCleavability.peptides_clev, input_peptides]),
        expr_col_name = expr_col_name,
        logtransform_expr = logtransform_expr,
        runtime_params = standard_runtime
    }
  }

  File final_peptide_list = select_first([ LogTransformExpr.peptides_expr,
                                           AggregatePeptides.peptides_agg,
                                           PredictCleavability.peptides_clev,
                                           input_peptides])
  ### Split merged peptides by length
  call SplitPeptidesPerLength {
    input:
      analysis_name = analysis_name,
      peptide_list = final_peptide_list,
      peptide_col_name = peptide_col_name,
      runtime_params = standard_runtime
  }

  ### Predict with NetMHC
  if (defined(netmhc_tar_gz)) {
    scatter (peps_allele in cross(SplitPeptidesPerLength.pep_len_files, ParseAlleles.alleles)) {
      call PredictNetMHC {
        input:
          netmhc_tar_gz = select_first([netmhc_tar_gz]),
          peptide_list = peps_allele.left,
          allele = peps_allele.right,
          peptide_col_name = peptide_col_name,
          runtime_params = standard_runtime
      }
    }
  }

  ### Generate peptide sequences features: dummy encoding only 
  ### (blosum and fuzzy generated on demand from the dummy encoding upon prediction)
  scatter (pepfilelen in SplitPeptidesPerLength.pep_len_files) {
    call FeaturizeEncoding {
      input:
        peptide_list = pepfilelen,
        encoding = "dummy",
        peptide_col_name = peptide_col_name,
        runtime_params = standard_runtime
    }
  }

  ### Predict with MS models (now includes ranks)
  scatter (allele_featfile in cross(ParseAlleles.alleles, FeaturizeEncoding.feature_files)) {
    call PredictMsAllele {
      input:
        patient_id = analysis_name,
        featfile = allele_featfile.right,
        exists_ctex = exists_ctex,
        exists_expr = exists_expr,
        run_netmhc = run_netmhc,
        allele = allele_featfile.left,
        models_path = models_path,
        peptide_col_name = peptide_col_name,
        features_all = features_all,
        feature_sets = feature_sets,
        runtime_params = standard_runtime
    }
  }

  Array[File] mspred_files = flatten(PredictMsAllele.mspred_files)

  ### Merge predictions into ensemble model scores: MSIntrinsic, MSIntrinsicC, MSIntrinsicEC; optinally merge NetMHC scores
  scatter (len in SplitPeptidesPerLength.pep_lens) {
    call MergeMspreds {
      input:
        patient_id = analysis_name,
        alleles = ParseAlleles.alleles,
        len = len,
        exists_ctex = exists_ctex,
        exists_expr = exists_expr,
        mspred_files = mspred_files,
        peptide_col_name = peptide_col_name,
        run_netmhc = run_netmhc,
        netmhc_EL = PredictNetMHC.netmhc_EL,
        netmhc_BA = PredictNetMHC.netmhc_BA,
        runtime_params = standard_runtime
    }
  }

  ### Convert scores to ranks
  scatter (preds_file in MergeMspreds.mspreds_file_wide) {
    call GetRanks {
      input:
        patient_id = analysis_name,
        alleles = ParseAlleles.alleles,
        lens = SplitPeptidesPerLength.pep_lens,
        exists_ctex = exists_ctex,
        exists_expr = exists_expr,
        run_netmhc = run_netmhc,
        models_path = models_path,
        peptide_col_name = peptide_col_name,
        preds_file = preds_file,
        runtime_params = standard_runtime
    }
  }

  ### Concatenate ranks files for all lengths
  call RanksConcatLengths {
    input:
      patient_id = analysis_name,
      peptide_col_name = peptide_col_name,
      exists_ctex = exists_ctex,
      exists_expr = exists_expr,
      run_netmhc = run_netmhc,
      assign_by_ranks_or_scores = assign_by_ranks_or_scores,
      assign_threshold = assign_threshold,
      assign_colors = assign_colors,
      ranks_len_files = GetRanks.ranks_file,
      wild_type_peptides_tar_gz = wild_type_peptides_tar_gz,
      runtime_params = standard_runtime
  }

  ### Workflow level outputs
  output {
    File all_binders = RanksConcatLengths.all_binders
    File mut_binders = RanksConcatLengths.mut_binders
    File allele_assignment_counts = RanksConcatLengths.allele_assignment_counts
    File allele_assignment_plots = RanksConcatLengths.allele_assignment_plots
  }
}


### Parse HLA alleles
task ParseAlleles {
  input {
    String analysis_name
    File alleles_file
    Runtime runtime_params
  }

  Int disk_gb = ceil(size(alleles_file,"GB")*2.5) + runtime_params.initial_disk_size

  command <<<
    set -e

    dos2unix ~{alleles_file}
    awk '{gsub("HLA-|*","")}1' ~{alleles_file} | \
    tr -d ":" | sort -u > parsed_alleles.txt

    python <<CODE
    import pandas as pd
    feature_file = '/data/feature_sets/ABCG_prot.parsed.clean.SEQ.ME.ALL.FEATS.txt'
    df = pd.read_csv(feature_file, sep=' ')
    alleles = set([x.rstrip() for x in open('parsed_alleles.txt')])
    # record HLA alleles that do not exist in the feature file
    fo = open('~{analysis_name}.removed_alleles.txt','w')
    alleles_not_exist = alleles.difference(df['allele'])
    fo.write('\n'.join(sorted(alleles_not_exist)))
    fo.close()
    # write only HLA alleles that exist in the feature file
    fo = open('~{analysis_name}.parsed_alleles.txt','w')
    fo.write('\n'.join(sorted(alleles.difference(alleles_not_exist))))
    fo.close()
    CODE

    { grep -c ^[A-C] ~{analysis_name}.removed_alleles.txt > removed.txt || true; }
  >>>

  output {
    Array[String] alleles = read_lines("~{analysis_name}.parsed_alleles.txt")
    File removed_alleles_file = "~{analysis_name}.removed_alleles.txt"
    Int removed_alleles = read_int("removed.txt")
  }
  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }
}


### Merge peptide lists if second list is supplied
task MergePeptideLists {
  input {
    String analysis_name
    File peptide_list
    File second_peptide_list
    Runtime runtime_params
  }

  Int disk_gb = ceil(size(peptide_list,"GB")+size(second_peptide_list,"GB"))*2 + runtime_params.initial_disk_size

  command <<<
    set -e

    dos2unix ~{peptide_list}
    cat ~{peptide_list} > ~{analysis_name}.merged_peps.txt
    tail -n+2 ~{second_peptide_list} >> ~{analysis_name}.merged_peps.txt
  >>>

  output {
    File merged_peptide_list = "~{analysis_name}.merged_peps.txt"
  }
  
  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }
}


### Run cleavability predictor
task PredictCleavability {
  input {
    String analysis_name
    File peptide_list
    String peptide_col_name
    Runtime runtime_params
  }

  Int disk_gb = ceil(size(peptide_list,"GB")*2) + runtime_params.initial_disk_size

  command <<<
    set -e

    dos2unix ~{peptide_list}

    Rscript /clevnn/featurize/featurize_clevnn.R \
      --infile ~{peptide_list} \
      --peptide_col_name ~{peptide_col_name} \
      --output clevnnfeats.tmp

    python /clevnn/pred/clevnn_pred.py \
      --input_file clevnnfeats.tmp \
      --output_file ~{analysis_name}.pep_clev.txt \
      --model_file /clevnn/pred/model_double20_drop0.2_features_375.txt_e10_bs10000_nh150.h5 \
      --features_file /clevnn/pred/features_375.txt
  >>>

  output {
    File peptides_clev = "~{analysis_name}.pep_clev.txt"
  }
    
  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }
}


### Aggregate per peptide
task AggregatePeptides {
  input {
    String analysis_name
    File peptide_list
    String peptide_col_name
    Boolean exists_expr
    String expr_col_name
    Runtime runtime_params
  }

  Int disk_gb = ceil(size(peptide_list,"GB")*2) + runtime_params.initial_disk_size

  command <<<
    set -e

    dos2unix ~{peptide_list}

    Rscript /aggregate/aggregate_pep.R \
      --peptide_col_name ~{peptide_col_name} \
      --exists_expr ~{exists_expr} \
      --expr_col_name ~{expr_col_name} \
      --infile ~{peptide_list} \
      --outfile ~{analysis_name}.pep_agg.txt
  >>>

  output {
    File peptides_agg = "~{analysis_name}.pep_agg.txt"
  }
  
  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }
}


### Ensure expression column is named log2TPM and optionally log-transform
task LogTransformExpr {
  input {
    String analysis_name
    File peptide_list
    String expr_col_name
    Boolean logtransform_expr
    Runtime runtime_params
  }

  Int disk_gb = ceil(size(peptide_list,"GB")*2) + runtime_params.initial_disk_size

  command <<<
    set -e

    dos2unix ~{peptide_list}

    Rscript /logtransform_expr/logtransform_expr.R \
      --infile ~{peptide_list} \
      --outfile ~{analysis_name}.pep_expr.txt \
      --expr_col_name ~{expr_col_name} \
      --logtransform_expr "~{logtransform_expr}"
  >>>

  output {
    File peptides_expr = "~{analysis_name}.pep_expr.txt"
  }
  
  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }
}


### Split merged peptides by length
task SplitPeptidesPerLength {
  input {
    String analysis_name
    File peptide_list
    String peptide_col_name
    Runtime runtime_params
  }

  Int disk_gb = ceil(size(peptide_list,"GB")*2.5) + runtime_params.initial_disk_size

  command <<<
    set -e

    echo peptide_list=~{peptide_list}
    dos2unix ~{peptide_list}

    python <<CODE
    peptide_fi = open('~{peptide_list}','r')
    header = peptide_fi.readline().rstrip('\n').split('\t')
    pep_col = -1
    try:
        pep_col = header.index('~{peptide_col_name}')
    except ValueError:
        raise Exception('~{peptide_col_name} does not exist in ~{peptide_list}')
    fo = dict()
    for line in peptide_fi:
        line = line.rstrip('\n').split('\t')
        pep_len = len(line[pep_col])
        if pep_len < 8 or 12 < pep_len:
            continue
        if pep_len not in fo:
            fo[pep_len] = open('~{analysis_name}.peps.{}.txt'.format(pep_len),'w')
            fo[pep_len].write('\t'.join(header)+'\n')
        fo[pep_len].write('\t'.join(line)+'\n')
    if fo:
        lens_fo = open('lens_present.txt','w')
        for pep_len in fo:
            lens_fo.write(str(pep_len)+'\n')
            fo[pep_len].close()
        lens_fo.close()
    else:
        raise Exception('~{peptide_list} does not contain 8, 9, 10, or 11aa peptides')
    CODE
  >>>

  output {
    Array[File] pep_len_files = glob("~{analysis_name}.peps.*.txt")
    Array[String] pep_lens = read_lines("lens_present.txt")
  }
  
  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }
}


### Predict with NetMHC
task PredictNetMHC {
  input {
    File netmhc_tar_gz
    File peptide_list
    String allele
    String peptide_col_name
    Runtime runtime_params
  }

  String prefix = basename(peptide_list)
  Int disk_gb = ceil(size(peptide_list,"GB")*2 + size(netmhc_tar_gz,"GB")*5) + runtime_params.initial_disk_size

  command <<<
    set -e

    ### Untar NetMHC executables and create a tmpdir
    NETMHC_DIR="/opt/netMHCpan"
    mkdir $NETMHC_DIR
    tar zxvf ~{netmhc_tar_gz} -C $NETMHC_DIR --strip-components 1
    mkdir -p $NETMHC_DIR/tmp

    allelePan="HLA-"$(echo ~{allele} | cut -c1-3)":"$(echo ~{allele} | cut -c4-5)
    echo NetMHC processing HLA allele: "$allelePan"

    ### NetMHC takes as input a simple file with only the peptides to predict for (skip the header line)
    pep_col_id=`head -1 ~{peptide_list} | sed -n $'1s/\t/\\\n/gp' | grep -nx ~{peptide_col_name} | cut -d: -f1`
    echo pep_col_id=$pep_col_id
    tail -n+2 ~{peptide_list} | awk -F'\t' -v i="$pep_col_id" '{print $i}' > ~{prefix}.tmp
    echo peponlyfile=~{prefix}.tmp

    ### Check if NetMHC supports the allele (the || true part prevents exit due to error due ot the -e)
    nmhc_support=`$NETMHC_DIR/Linux_x86_64/bin/netMHCpan -listMHC | grep ~{allele} || true`
    echo nmhc_support = $nmhc_support
    if [ -z "$nmhc_support" ]; then
      echo  "Warning: NetMHC does not support allele: ~{allele}"
        ### Create pseudo-netmhc output files with scores and ranks set to -1 / 100 accordingly
        ### (this accomodates a mix of supported and unsupported alleles)
        ### Run with a fixed allele, then edit columns accordingly
        allelePan_dummy="HLA-A01:01"
        allelePan=$allelePan_dummy
    fi

    echo allelePan: $allelePan

    ### EL
    echo Running netMHCPan EL
    $NETMHC_DIR/Linux_x86_64/bin/netMHCpan -a "$allelePan" -p ~{prefix}.tmp -xls -xlsfile \
      ~{prefix}.netmhcpan40.~{allele}.EL.xls.txt > ~{prefix}.netmhcpan40.~{allele}.EL

    ### BA - note -BA is for binding affinity prediction (i.e. nM, while the default is eluted ligand likelihood EL)
    echo Running netMHCPan BA
    $NETMHC_DIR/Linux_x86_64/bin/netMHCpan -BA -a "$allelePan" -p ~{prefix}.tmp -xls -xlsfile \
      ~{prefix}.netmhcpan40.~{allele}.BA.xls.txt > ~{prefix}.netmhcpan40.~{allele}.BA

    ### Trim footer (last 5 lines) from output files with ranks are not truncated to 4 digits
    head -n -5 ~{prefix}.netmhcpan40.~{allele}.EL > ~{prefix}.netmhcpan40.~{allele}.EL.trim
    head -n -5 ~{prefix}.netmhcpan40.~{allele}.BA > ~{prefix}.netmhcpan40.~{allele}.BA.trim

    ### Check NetMHCPan version (it could be NetMHCPan-4.1 although the suffix is 'netmhcpan40')
    version41=`grep "NetMHCpan version 4.1" ~{prefix}.netmhcpan40.~{allele}.BA.trim | wc -l`
    if [ $version41 == 1 ]; then
        ### Remove EL columns from BA
        awk '{$12=""; $13=""; print $0}' ~{prefix}.netmhcpan40.~{allele}.BA.trim > ~{prefix}.netmhcpan40.~{allele}.BA.trim.tmp
        mv ~{prefix}.netmhcpan40.~{allele}.BA.trim.tmp ~{prefix}.netmhcpan40.~{allele}.BA.trim
    
        ### Extract headers
        awk 'BEGIN{ hline=0 }{ if(substr($1,1,3)=="---"){hline+=1} print $0; if(hline==2){exit 0} }' ~{prefix}.netmhcpan40.~{allele}.EL.trim > ~{prefix}.netmhcpan40.~{allele}.EL.trim.tmp
        awk 'BEGIN{ hline=0 }{ if(substr($1,1,3)=="---"){hline+=1} print $0; if(hline==2){exit 0} }' ~{prefix}.netmhcpan40.~{allele}.BA.trim > ~{prefix}.netmhcpan40.~{allele}.BA.trim.tmp
        
        ### Append lines to make the number of lines same as NetMHCPan-4.0
        lc=`cat ~{prefix}.netmhcpan40.~{allele}.BA.trim.tmp | wc -l`
        echo '' | awk -v h=$lc 'BEGIN{for(i=1; i<=51-h; ++i){print "#"} }' | cat - ~{prefix}.netmhcpan40.~{allele}.BA.trim > ~{prefix}.netmhcpan40.~{allele}.BA.trim.tmp
        lc=`cat ~{prefix}.netmhcpan40.~{allele}.EL.trim.tmp | wc -l`
        echo '' | awk -v h=$lc 'BEGIN{for(i=1; i<=50-h; ++i){print "#"} }' | cat - ~{prefix}.netmhcpan40.~{allele}.EL.trim > ~{prefix}.netmhcpan40.~{allele}.EL.trim.tmp
        
        ### Rename formatted files
        mv ~{prefix}.netmhcpan40.~{allele}.BA.trim.tmp ~{prefix}.netmhcpan40.~{allele}.BA.trim
        mv ~{prefix}.netmhcpan40.~{allele}.EL.trim.tmp ~{prefix}.netmhcpan40.~{allele}.EL.trim
    fi

    ### Cleanup
    rm ~{prefix}.netmhcpan40.~{allele}.EL
    rm ~{prefix}.netmhcpan40.~{allele}.BA

    if [ -z "$nmhc_support" ]; then    
      ### Edit columns
      sed -i.bak -e "s/HLA-A[*]\?01:01/HLA_dummy/g" ~{prefix}.netmhcpan40.~{allele}.EL.trim
      sed -i.bak -e "s/HLA-A[*]\?01:01/HLA_dummy/g" ~{prefix}.netmhcpan40.~{allele}.BA.trim
      head -n 50 ~{prefix}.netmhcpan40.~{allele}.EL.trim > ~{prefix}.netmhcpan40.~{allele}.EL.trim.tmp
      head -n 51 ~{prefix}.netmhcpan40.~{allele}.BA.trim > ~{prefix}.netmhcpan40.~{allele}.BA.trim.tmp
      tail -n +51 ~{prefix}.netmhcpan40.~{allele}.EL.trim | \
        awk '{$12=0; $13=100; print}' >> ~{prefix}.netmhcpan40.~{allele}.EL.trim.tmp
      tail -n +52 ~{prefix}.netmhcpan40.~{allele}.BA.trim | \
        awk '{$12=0; $13=50000; $14=100; print}' >> ~{prefix}.netmhcpan40.~{allele}.BA.trim.tmp
      mv ~{prefix}.netmhcpan40.~{allele}.EL.trim.tmp ~{prefix}.netmhcpan40.~{allele}.EL.trim
      mv ~{prefix}.netmhcpan40.~{allele}.BA.trim.tmp ~{prefix}.netmhcpan40.~{allele}.BA.trim

      ### Cleanup
      rm ~{prefix}.netmhcpan40.~{allele}.EL.trim.bak
      rm ~{prefix}.netmhcpan40.~{allele}.BA.trim.bak
    fi

    ### Cleanup
    rm ~{prefix}.tmp
  >>>

  output {
    # keep these (over the xlsx files) because ranks are not truncated to 4 digits
    File netmhc_EL = "~{prefix}.netmhcpan40.~{allele}.EL.trim"
    File netmhc_BA = "~{prefix}.netmhcpan40.~{allele}.BA.trim"
  }
  
  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }
}


### Generate peptide sequences features: dummy, blosum, or fuzzy encoding
task FeaturizeEncoding {
  input {
    File peptide_list
    String encoding
    String peptide_col_name
    Runtime runtime_params
  }

  String prefix = basename(peptide_list, ".txt")

  Int disk_gb = ceil(size(peptide_list,"GB")*2.5) + runtime_params.initial_disk_size

  command <<<
    set -e
    
    Rscript /encoding/encoding.R \
      --infile ~{peptide_list} \
      --outfile ~{prefix}.~{encoding} \
      --peptide_col_name ~{peptide_col_name} \
      --coding ~{encoding}
  >>>

  output {
    File feature_files = "~{prefix}.~{encoding}"
  }
  
  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }
}


### Predict with MS models (either allele-and-len-specific model, or panpan model when missing)
task PredictMsAllele {
  input {
    String patient_id
    File featfile
    Boolean exists_ctex
    Boolean exists_expr
    Boolean run_netmhc
    String allele
    String models_path
    String peptide_col_name
    String features_all
    String feature_sets
    Runtime runtime_params
  }

  Int disk_gb = ceil(size(featfile,"GB")*3) + runtime_params.initial_disk_size

  command <<<
    set -e

    export TMPDIR=/tmp # to fix error: AF_UNIX path too long

    ### Get the len from the input file name    
    len=`basename ~{featfile} | tr "." "\n" | grep "pep" -A1 | tail -n1`
    echo len=$len

    ### Figure out whether to use allele-and-length specific or panpan model
    ### ugly shenanigans splitting grep in two to make wdl happy
    echo grep -w ~{allele} /data/model_exists_dat.txt | grep -w $len | cut -d ' ' -f 3
    use_specific=`grep -w ~{allele} /data/model_exists_dat.txt | grep -w $len | cut -d ' ' -f 3`
    echo use_specific = $use_specific

    ### Copy models from google bucket and extract the tar.gz file
    if [ "$use_specific" == 1 ]; then
      echo using allele-and-length specific model
      local_model_dir=./~{allele} # created by the untar
      gsutil cp "~{models_path}~{allele}.tar.gz" .
      tar -xzf ~{allele}.tar.gz .

      models_linear_RDS_file=models_linear_"$len"_~{allele}_slim.RDS
      gsutil cp ~{models_path}models_linear/$models_linear_RDS_file "$local_model_dir/$models_linear_RDS_file"
    else
      echo using panpan model
      local_model_dir="./panpan"
      mkdir $local_model_dir # in the specific case, the untar creates the allele-specific model folder
      gsutil -m cp -r ~{models_path}models_panpan/models_pan_pan_CV/ "$local_model_dir"

      hla=`echo ~{allele} | cut -c1`
      models_linear_RDS_file=models_linear_"$len"_"$hla"_slim.RDS
      gsutil cp ~{models_path}models_panpan/models_linear_pan_pan/"$models_linear_RDS_file" "$local_model_dir/$models_linear_RDS_file"
    fi

    ### Predict MSi (all encodings at once, also only one feature_set since the rest are linear models on top)
    feature_set=~{feature_sets}
    if [ "$use_specific" == 1 ]; then
      ### Execute allele-and-len-specific model
      python /predict/mspredict/ann_pred_py3_optimize.py --sample_name ~{patient_id} --alleles ~{allele} --len $len \
        --peptide_col_name ~{peptide_col_name} --features_pep_all ~{features_all} --feature_set_pep $feature_set \
        --model_path "./" --infile ~{featfile} --outpath "./"
    else
      ### Execute panpan model
      python /predict/mspredict/ann_pred_py3_optimize_pan.py --sample_name ~{patient_id} --alleles ~{allele} --len $len \
        --peptide_col_name ~{peptide_col_name} --features_pep_all ~{features_all} --feature_set_pep $feature_set \
        --feature_set_hla_all "allele_feats" --feature_set_hla "allele_feats_PCA" \
        --model_path "./panpan" --infile ~{featfile} --outpath "./"
    fi

    ### Predict MSiC, MSiCE, MSiCEB (linear models)
    if [ "~{exists_ctex}" = true ] || [ "~{exists_expr}" = true ]; then
      infile=./mspred_"~{patient_id}"_"~{allele}"_"$len"_"~{feature_sets}"_summary.txt
      outfile=./mspred_"~{patient_id}"_"~{allele}"_"$len"_"~{feature_sets}"_summary_L.txt
      Rscript /predict/msevallinear/mseval_linear.R \
        --len $len --allele ~{allele} \
        --models_linear_RDS_file $local_model_dir/$models_linear_RDS_file \
        --exists_ctex ~{exists_ctex} --exists_expr ~{exists_expr} \
        --infile $infile --outfile $outfile
      mv $outfile $infile
    fi

    ### Cleanup
    if [ -f ./~{allele}.tar.gz ]; then
      rm -r ~{allele}.tar.gz ./~{allele}
    else
      rm -r ./panpan
    fi
  >>>

  output {
    Array[File] mspred_files = glob('*_summary.txt')
  }
  
  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }
}


### Merge MS predictions optionally with NetMHC scores
task MergeMspreds {
  input {
    String patient_id
    Array[String] alleles
    String len
    Boolean exists_ctex
    Boolean exists_expr
    String peptide_col_name
    Array[File] mspred_files
    Boolean run_netmhc
    Array[File]? netmhc_EL
    Array[File]? netmhc_BA
    
    Runtime runtime_params
  }

  Int disk_gb = ceil(size(mspred_files,"GB")*3) + runtime_params.initial_disk_size

  command <<<
    ### Get the index of the peptide column
    echo peptide_col_name=~{peptide_col_name}
    echo len=~{len}
    output_file_name="mspreds.~{patient_id}.pep.~{len}"
    echo output_file_name=$output_file_name

    Rscript /predict/msmerge/msmerge_optimized.R --patient ~{patient_id} \
      --alleles ~{sep=',' alleles} --len ~{len} --exists_ctex ~{exists_ctex} \
      --exists_expr ~{exists_expr} --peptide_col_name ~{peptide_col_name} \
      --infiles ~{sep=',' mspred_files} --outfile $output_file_name \
      --run_netmhc ~{run_netmhc} --netmhc_EL "~{sep=',' netmhc_EL}" \
      --netmhc_BA "~{sep=',' netmhc_BA}"
  >>>

  output {
    File mspreds_file_long = "mspreds.~{patient_id}.pep.~{len}.long.txt"
    File mspreds_file_wide = "mspreds.~{patient_id}.pep.~{len}.wide.txt"
  }
  
  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }
}


task GetRanks {
  input {
    String patient_id
    Array[String] alleles
    Array[String] lens
    Boolean exists_ctex
    Boolean exists_expr
    Boolean run_netmhc
    String models_path
    String peptide_col_name
    File preds_file
    Runtime runtime_params
  }

  Int disk_gb = ceil(size(preds_file,"GB")*2.5) + runtime_params.initial_disk_size

  command <<<
    set -e

    export TMPDIR=/tmp # to fix error: AF_UNIX path too long
    len=`basename ~{preds_file} | tr "." "\n" | grep "pep" -A1 | tail -n1`
    echo len=$len

    ### Copy ecdf files from google bucket
    for allele in ~{sep=' ' alleles}; do
      ### Figure out whether to use allele-and-length specific or panpan model
      ### ugly shenanigans splitting grep in two to make wdl happy
      echo grep -w "$allele" /data/model_exists_dat.txt | grep -w $len | cut -d ' ' -f 3
      use_specific=`grep -w "$allele" /data/model_exists_dat.txt | grep -w $len | cut -d ' ' -f 3`
      echo use_specific = $use_specific

      if [ "$use_specific" == 1 ]; then
        echo "$allele": using allele-and-length specific model
        model_ecdf=ecdf_panpan_"$allele"_nmhc_mhcflurry.RDS # 'panpan' here is just a typo, these are specific ecdfs
        gsutil cp ~{models_path}ecdf/"$model_ecdf" .
      else
        echo "$allele": using panpan model
        # most ecdfs only have MSiX scores, some (the 95+some) have other tools as well that were used in the eval
        model_ecdf=ecdf_panpan_"$allele"_nmhc_mhcflurry.RDS
        set +e # turn off option temporarily because despite the -q an error code is returned if file not found!
        gsutil -q stat ~{models_path}models_panpan/ecdf/"$model_ecdf"
        return_value=$? # Result 0 means exists; 1 means not exists
        set -e # turn back on
        if [ $return_value == 1 ]; then
          echo ~{models_path}models_panpan/ecdf/"$model_ecdf" not found, using ecdf_panpan_"$allele".RDS
          model_ecdf=ecdf_panpan_"$allele".RDS
        fi
        gsutil cp ~{models_path}models_panpan/ecdf/"$model_ecdf" .
      fi
    done

    ### Ranks
    Rscript /predict/ecdf/ecdf_wide.R \
      --len $len --alleles ~{sep=',' alleles} \
      --exists_ctex ~{exists_ctex} --exists_expr ~{exists_expr} --run_netmhc ~{run_netmhc} \
      --ecdf_model_path "./" \
      --infile ~{preds_file} --outfile ~{patient_id}.ranks
  >>>

  output {
    File ranks_file = "~{patient_id}.ranks"
  }
  
  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }
}


### Concatenate ranks files for all lengths, assign alleles, provide summaries
task RanksConcatLengths {
  input {
    String patient_id
    String peptide_col_name
    Boolean exists_ctex
    Boolean exists_expr
    Boolean run_netmhc
    String assign_by_ranks_or_scores
    Float assign_threshold
    String assign_colors
    Array[File] ranks_len_files
    File? wild_type_peptides_tar_gz
    Runtime runtime_params
  }
  
  String outfile = "~{patient_id}.ranks.cat"
  Int disk_gb = ceil(size(ranks_len_files,"GB")*2.5) + runtime_params.initial_disk_size

  ### using <<< syntax due to parsing problems with {} and awk expression
  command <<<
    set -e

    # Only keep header from first file:
    #    FNR is the number of lines (records) read so far in the current file.
    #     NR is the number of lines read overall.
    awk 'FNR==1 && NR!=1{next;}{print}' ~{sep=" " ranks_len_files} > ~{outfile}

    # Assign peptides to alleles and make summary plots
    Rscript /predict/assign/assign.R \
      --sample_name ~{patient_id} --peptide_col_name ~{peptide_col_name} \
      --exists_ctex ~{exists_ctex} --exists_expr ~{exists_expr} --run_netmhc ~{run_netmhc} \
      --assign_by_ranks_or_scores ~{assign_by_ranks_or_scores} --assign_threshold ~{assign_threshold} \
      --assign_colors "~{assign_colors}" --infile ~{outfile} --outpath "./"

    # append NeoORF information
    python3 /root/neovax-pipeline/add_neoorf_status.py ~{patient_id}_predictions.txt ~{patient_id}_predictions.neoorf.txt

    if [[ -f "~{wild_type_peptides_tar_gz}" ]] ; then
      mkdir wild_type_pep_dir
      tar zxvf ~{wild_type_peptides_tar_gz} -C wild_type_pep_dir --strip-components 1

      # filter wild type peptide
      python3 /root/neovax-pipeline/filter_peptide.py wild_type_pep_dir ~{patient_id}_predictions.neoorf.txt
      mv ~{patient_id}_predictions.neoorf.txt.only_non_wild_type.tsv ~{patient_id}.mut_binders.txt
      mv ~{patient_id}_predictions.neoorf.txt.wild_type_status.tsv ~{patient_id}.all_binders.txt
    else
      mv ~{patient_id}_predictions.neoorf.txt ~{patient_id}.all_binders.txt
      touch ~{patient_id}.mut_binders.txt
    fi
    mv ~{patient_id}_allele_assignment_counts.txt ~{patient_id}.allele_assignment_counts.txt
    mv ~{patient_id}.pdf ~{patient_id}.allele_assignment_plots.pdf
  >>>

  output {
    File all_binders = "~{patient_id}.all_binders.txt"
    File mut_binders = "~{patient_id}.mut_binders.txt"
    File allele_assignment_counts = "~{patient_id}.allele_assignment_counts.txt"
    File allele_assignment_plots = "~{patient_id}.allele_assignment_plots.pdf"
  }
  
  runtime {
    docker        : runtime_params.docker
    bootDiskSizeGb: runtime_params.boot_disk_size
    preemptible   : runtime_params.preemptible
    cpu           : runtime_params.cpu
    disks         : "local-disk " + disk_gb + " HDD"
    memory        : runtime_params.mem + "GB"
    maxRetries    : runtime_params.max_retries
  }
}