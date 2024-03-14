version 1.0

workflow FormatHLA {
  input {
    String patient_id
    File optitype_hla
    File polysolver_hla
    File? neovax_clinical_hla
  }

  call FormatPredictedHLA {
    input:
      patient_id = patient_id,
      optitype_hla = optitype_hla,
      polysolver_hla = polysolver_hla
  }
  
  if(defined(neovax_clinical_hla)) {
    call FormatClinicalHLA {
      input:
        patient_id = patient_id,
        clinical_hla = select_first([neovax_clinical_hla])
    }
  }

  output {
    File? clinical_hla = FormatClinicalHLA.formatted_hla
    File predicted_hla = FormatPredictedHLA.formatted_hla
    File predicted_hla_status = FormatPredictedHLA.concordance
  }
}


task FormatPredictedHLA {
  input {
    String patient_id
    File optitype_hla
    File polysolver_hla
  }

  command <<<
    python <<CODE
    # Load formatted HLA alleles
    optitype = [x.rstrip() for x in open('~{optitype_hla}')]
    polysolver = [x.rstrip() for x in open('~{polysolver_hla}')]
    # Stop if there are no readable alleles
    if len(optitype) == 0 and len(polysolver) == 0:
        raise Exception('No readable alleles!')
    fo = open('~{patient_id}.prediction_status.txt', 'w')
    if sorted(optitype) == sorted(polysolver):
        fo.write('concordant\n')
    else:
        fo.write('discordant\n')
    fo.close()
    # Record non-redundant HLA alleles (but keep homozygous alleles)
    predicted = optitype[:]
    for x in polysolver:
        if x not in predicted:
            predicted.append(x)
    fo = open('~{patient_id}.predicted_hla.txt', 'w')
    fo.write('\n'.join(sorted(predicted)))
    fo.close()
    CODE
  >>>

  output {
    File formatted_hla = "~{patient_id}.predicted_hla.txt"
    String concordance = read_string("~{patient_id}.prediction_status.txt")
  }

  runtime {
    docker: "python:3"
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    preemptible: 1
    maxRetries: 1
  }
  
  meta {
    author: "Junko Tsuji"
  }
}


task FormatClinicalHLA {
  input {
    String patient_id
    File clinical_hla
  }
  
  command <<<
    python <<CODE
    import re
    import itertools
    allele_counts = 0
    fo = open('~{patient_id}.clinical_hla.txt', 'w')
    hla_alleles = [h.replace('*','').split() for h in open('~{clinical_hla}')]
    hla_alleles = list(itertools.chain(*hla_alleles))
    for allele in hla_alleles:
        # Write original line if already formatted
        if re.match('^HLA-[ABC][0-9]{2}:[0-9]{2}$', allele):
            fo.write(allele+'\n')
            allele_counts += 1
        elif allele[0] in ['A', 'B', 'C']:
            fo.write('HLA-'+allele+'\n')
            allele_counts += 1
    fo.close()
    # Stop if there are no readable alleles
    if allele_counts == 0:
        raise Exception('No readable alleles!')
    CODE
  >>>

  output {
    File formatted_hla = "~{patient_id}.clinical_hla.txt"
  }

  runtime {
    docker: "python:3"
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    preemptible: 1
    maxRetries: 1
  }
  
  meta {
    author: "Junko Tsuji"
  }
}
