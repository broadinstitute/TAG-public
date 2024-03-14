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

workflow LegoPlotter {
  input {
    String pair_name
    File maf_file
    File mut_categs
    
    # 'exome', 'genome', or 'unit' to pick trimer coverage for normalization
    # 'exome' = typical exome coverage
    # 'genome' = genome coverage
    # 'unit' = equal weight across bases
    String cov_string = "exome"

    # workflow level runtime parameters
    Int preemptible = 2
    Int max_retries = 1
    Int additional_disk = 20
    Int boot_disk_size = 15
    Int mem = 8
    Int cpu = 1
  }

  Runtime standard_runtime = { "preemptible": preemptible,
                               "max_retries": max_retries,
                               "mem": mem,
                               "cpu": cpu,
                               "docker": "us.gcr.io/tag-public/neovax-tag-cga:v1",
                               "boot_disk_size": boot_disk_size,
                               "initial_disk_size": additional_disk }

	call CreateLegoPlot {
    input:
      pair_name = pair_name,
      maf_file = maf_file,
      mut_categs = mut_categs,
      cov_string = cov_string,
      runtime_params = standard_runtime
  }
  
  output {
    Array[File] lego_plotter_ais = CreateLegoPlot.ais
    Array[File] lego_plotter_pngs = CreateLegoPlot.pngs
    Array[File] lego_plotter_figs = CreateLegoPlot.figs
    Array[File] lego_plotter_pss = CreateLegoPlot.pss
    File mut_legos_html = CreateLegoPlot.mut_legos_html
  }
}

task CreateLegoPlot {
  input {
    String pair_name
    File maf_file
    File mut_categs
    String cov_string
    Runtime runtime_params
  }
  
  Int disk_gb = runtime_params.initial_disk_size + ceil(size(maf_file,"GB") + size(mut_categs, "GB"))

  command <<<
    set -e

    # create directory where category data is expected
    mkdir -pv /xchip/cga/reference/hg19/annotation/db/tracks/hg19/c65e/

    # copy the category data
    cp -vf ~{mut_categs} /xchip/cga/reference/hg19/annotation/db/tracks/hg19/c65e/

    # correct [Start|End]_Position to [Start|End]_position
    python /usr/local/bin/subset_maf_columns.py --correct-maf-dialect \
      /usr/local/bin/maf_columns/basic_maf_columns.txt ~{maf_file} ~{pair_name}.corrected.maf

    # split MNPs/ONPs into SNPs
    python /usr/local/bin/split_onp_to_snp.py ~{pair_name}.corrected.maf ~{pair_name}.maf

    # run the plotter and acquire the exit code (plot rates and zscale)
    /usr/local/bin/run_call_lego_plotter.sh ~{pair_name}.maf . ~{cov_string} yes yes
    PLOTTER_EXIT_CODE=$? ;

    # encode the PNGS
    for PNG in `find *.png`; do 
      uuencode  -m $PNG  /dev/stdout | grep -Pv '^begin'|grep -Pv '^====$'|tr -d "\n" > $PNG.enc ; 
      UUE=$? ; 
      if [ "$UUE" -ne "0" ] ; then exit 2 ; fi ;
    done 

    # update the HTML in memory/place
    python <<CODE
    import glob
    import re
    for H in glob.glob('*_legos.html'):
        reader=open(H,'r')
        # store the HTML into a single string
        all_lines=''
        for line in reader:
            all_lines=all_lines+line
        reader.close()
        #iterate over the encs
        for ENC in glob.glob('*.enc'):
            print 'got enc = ',ENC
            PNG=re.sub('\.enc$','',ENC)
            enc_reader=open(ENC,'r')
            enc_data=''
            for enc_line in enc_reader:
                enc_data=enc_data+enc_line.strip()
            enc_reader.close()
            # embed the img, close the src, and add a title
            enc_replacement='data:image/png;base64,'+str(enc_data)+'\" title=\"'+str(PNG)+'\" '
            all_lines=all_lines.replace(PNG,enc_replacement)
        writer=open(H,'w')
        writer.write(all_lines)
        writer.close()
    CODE
    PY_EXIT=$? ;
    if [ "$PY_EXIT" -ne "0" ] ; then exit 3 ; fi ;

    # propagate exit code
    bash -c "exit $PLOTTER_EXIT_CODE"
  >>>

  output {
    Array[File] ais = glob("*.ai")
    Array[File] pngs = glob("*.png")
    Array[File] figs = glob("*.fig")
    Array[File] pss = glob("*.ps")
    File mut_legos_html = "~{pair_name}.maf.mutation_legos.html"
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

  meta {
    author: "Junko Tsuji"
    description: "Adopted from [the CGA WES Characterization Workflow](https://firecloud.terra.bio/#workspaces/broad-fc-getzlab-workflows/CGA_WES_Characterization_OpenAccess) and tweaked for handling ONPs and MNPs"
  }
}
