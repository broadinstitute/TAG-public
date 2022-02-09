workflow GenerateReportMMRF {
    String patient_id
    String sample_id
    String assay_version

    # Changing pf_bases_aligned to String from Int since
    # it appears Cromwell cannot handle a long int
    String pf_bases_aligned

    Float pct_selected_bases
    Int target_territory
    Int mean_duplex_depth
    File contamination_file

    Int pass_snv_number
    Int pass_indel_number

    # Signing an existing report
    File? input_report_html
    String? signature
    Boolean sign_and_date_report = defined(input_report_html) && defined(signature)

    if(sign_and_date_report) {
      call SignAndDateReport {
        input:
          out_basename = patient_id + "-" + sample_id,
          report_html = select_first([input_report_html]),
          signature = select_first([signature])
      }
    }
    if(!sign_and_date_report) {
      call ExtractContamination {
        input:
          contamination_table = contamination_file
      }
      call GenerateReportTemplate {
        input:
          patient_id = patient_id,
          sample_id = sample_id,
          assay_version = assay_version,
          pct_selected_bases = pct_selected_bases,
          target_territory = target_territory,
          pf_bases_aligned = pf_bases_aligned,
          mean_duplex_depth = mean_duplex_depth,
          contamination_file = contamination_file,
          pass_snv_number = pass_snv_number,
          pass_indel_number = pass_indel_number,
          out_basename = patient_id + "-" + sample_id
      }
    }
    call CompileReport {
      input:
        out_basename = patient_id + "-" + sample_id,
        report_html = select_first([SignAndDateReport.signed_report_html, GenerateReportTemplate.report_html])
    }
    output {
      Float? contam_frac = ExtractContamination.contam_frac
      File report_html = select_first([SignAndDateReport.signed_report_html, GenerateReportTemplate.report_html])
      File report_pdf = CompileReport.report_pdf
    }
}


task ExtractContamination {
   File contamination_table
   
   command <<<
      grep -v ^sample ${contamination_table} | awk '{print $2}' > contam.txt   
   >>>
   runtime {
      docker: "us.gcr.io/broad-gatk/gatk:4.1.9.0"
      memory: "2 GB"
      disk: "local-disk 30 HDD"
      maxRetries: 1
      preemptible: 3
   }
   output {
      Float contam_frac = read_float("contam.txt")
   }
}


task CompileReport {
    String out_basename
    File report_html

    command <<<
      set -e
      
      cp /usr/tag/template/broad_crsp_logo.png .

      wkhtmltopdf --page-size "Letter" --margin-bottom "10mm" \
                  --footer-html /usr/tag/template/report_footer.html --footer-line \
                  --header-html /usr/tag/template/report_header.html --enable-local-file-access \
                  "${report_html}" "${out_basename}-mmrf_technical_report.pdf"
    >>>
    output {
      File report_pdf = "${out_basename}-mmrf_technical_report.pdf"
    }
    runtime {
        docker: "us.gcr.io/tag-team-160914/mmrf-report-generation:v2"
        preemptible: 3
        maxRetries: 1
        disks: "local-disk 50 HDD"
        memory: "2 GB"    
    }
}


task SignAndDateReport {
    String out_basename
    File report_html
    String signature

    command <<<
      set -e

      # Set a signed date to the day the job is run
      DATE=`TZ=":US/eastern" date +"%m/%d/%Y"`

      PYTHON_CMD="
from datetime import date
fo = open('${out_basename}-mmrf_technical_report.html', 'w')
for line in open('${report_html}', 'r'):
    line = line.rstrip('\n')
    if line.startswith('<td style=\"border:none;width:100px;padding:0;border-bottom: 1pt solid black;\">'):
        line = line.split('</td>')[0] + ' <font size=\"2pt\"><i> ${signature} </i></font></td>'
    if line.startswith('<td style=\"border:none;width:60px;text-align:left;padding:0;border-bottom: 1pt solid black;\">'):
        line = line.split('</td>')[0] + ' <font size=\"2pt\"> $DATE </font></td>'
    fo.write(line + '\n')
fo.close()
"
      python3 -c "$PYTHON_CMD"
    >>>
    output {
        File signed_report_html = "${out_basename}-mmrf_technical_report.html"
    }
    runtime {
        docker: "us.gcr.io/tag-team-160914/mmrf-report-generation:v2"
        preemptible: 3
        maxRetries: 1
        disks: "local-disk 50 HDD"
        memory: "2 GB"    
    }
}


task GenerateReportTemplate {
    String out_basename
    String patient_id
    String sample_id
    String assay_version
    
    Float pct_selected_bases
    Int target_territory
    String pf_bases_aligned
    Int mean_duplex_depth
    File contamination_file

    Int pass_snv_number
    Int pass_indel_number

    command <<<
      set -e

      # Patient information
      echo "patient_id==${patient_id}" > input.txt
      echo "sample_id==${sample_id}" >> input.txt

      # Coverage information
      mean_raw_depth=`echo | awk -v b="${pf_bases_aligned}" \
                                 -v s="${pct_selected_bases}" \
                                 -v t="${target_territory}" '{print int( b*s/t )}'`
      echo "mean_raw_depth==$mean_raw_depth" >> input.txt
      echo "mean_duplex_depth==${mean_duplex_depth}" >> input.txt

      # PASS or FAIL information
      if [[ $mean_raw_depth -gt 57000 && ${mean_duplex_depth} -gt 400 ]];
      then
        echo "pass_or_fail==<span style='color:green;'>passed</span>" >> input.txt
      else
        echo "pass_or_fail==<span style='color:red;'>failed</span>" >> input.txt
      fi

      # Contamination information
      contamination=`tail -n1 ${contamination_file} | cut -f2 | awk '{printf "%.3f", $0*100}'`
      echo "contamination==$contamination" >> input.txt

      # Detected mutations
      echo "pass_snv_num==${pass_snv_number}" >> input.txt
      echo "pass_indel_num==${pass_indel_number}" >> input.txt

      # Assay version and time stamp
      echo "test_version==${assay_version}" >> input.txt
      TZ=":US/eastern" date +"%m/%d/%Y" | awk '{print "analysis_date=="$0}' >> input.txt

      # Run python code snippet to generate a sample specific HTML file
      PYTHON_CMD="
from string import Template
var_list = dict([tuple(i.rstrip().split('==')) for i in open('input.txt')])
template_html = Template(open('/usr/tag/template/report_template.html').read())
output_html = template_html.substitute(var_list)
    
fo = open('${out_basename}-mmrf_technical_report.html', 'w')
fo.write(output_html)
fo.close()
"
      python3 -c "$PYTHON_CMD"
    >>>
    output {
      File report_html = "${out_basename}-mmrf_technical_report.html"
    }
    runtime {
        docker: "us.gcr.io/tag-team-160914/mmrf-report-generation:v2"
        preemptible: 3
        maxRetries: 1
        disks: "local-disk 50 HDD"
        memory: "2 GB"
    }
}
