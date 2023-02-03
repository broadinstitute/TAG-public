task generateLegoPlot {
    File input_file
    String input_file_format
    String output_prefix
    String sample_id
    Boolean is_exome

    File? plotter_override
    File? renderer_override
    File? callable_regions
    File? ref_fasta

    String plotter_docker
    Int? mem
    Int? disk_space_gb
    Int? preepmtible_attempts

    String precomputed_option = if is_exome then "--mutsig-exome" else "--mutsig-genome"

    command <<<

        export PLOTTER_SRC=${default="/usr/tag/scripts/lego-plot.py" plotter_override}
        export RENDERER_SRC=${default="/usr/tag/scripts/lego-report.py" renderer_override}

        # Mutation rate spectrum
        python $PLOTTER_SRC --plot-title "${sample_id}: MutSig 2CV precomputed callable regions" \
                            --output-prefix precomputed_rate \
                            ${precomputed_option} \
                            ${"-s " + ref_fasta} \
                            ${input_file_format} ${input_file}
        if [[ -f "${callable_regions}"  ]]; then
           python $PLOTTER_SRC --plot-title "${sample_id}: Sample callable regions" \
                               --output-prefix sample_rate \
                               --user-coverage ${callable_regions} \
                               ${"-s " + ref_fasta} \
                               ${input_file_format} ${input_file}
        fi

        # Mutation count spectrum
        python $PLOTTER_SRC --plot-title "${sample_id}: All variants" \
                            --all-variants \
                            --output-prefix all_count \
                            ${"-s " + ref_fasta} \
                            ${input_file_format} ${input_file}
        python $PLOTTER_SRC --plot-title "${sample_id}: PASSed variants" \
                            --output-prefix pass_count \
                            ${"-s " + ref_fasta} \
                            ${input_file_format} ${input_file}
        
        # MAF ONLY: Mutation count spectrum sliced by allele fraction
        if [[  "${input_file_format}" == "maf" ]]; then
            python $PLOTTER_SRC --plot-title "0 <= AF < 0.1" \
                                --output-prefix af_0_01 \
                                --af-slice 0 0.1 maf ${input_file}
            python $PLOTTER_SRC --plot-title "0.1 <= AF < 0.25" \
                                --output-prefix af_01_025 \
                                --af-slice 0.1 0.25 maf ${input_file}
            python $PLOTTER_SRC --plot-title "0.25 <= AF < 0.5" \
                                --output-prefix af_025_05 \
                                --af-slice 0.25 0.5 maf ${input_file}
            python $PLOTTER_SRC --plot-title "0.5 <= AF < 1" \
                                --output-prefix af_05_1 \
                                --af-slice 0.5 1 maf ${input_file}
            ALLELE_SLICE_PDF="--allele-slice af_0_01.pdf af_025_05.pdf af_01_025.pdf af_05_1.pdf"
        fi

        # Summarize lego plots into slides
        python $RENDERER_SRC --output-prefix ${output_prefix} \
                             --mutation-rate `ls *_rate.pdf` \
                             --mutation-count all_count.pdf pass_count.pdf \
                             $ALLELE_SLICE_PDF
        pdflatex ${output_prefix}.tex
    >>>
    runtime {
        docker: plotter_docker
        memory: select_first([mem, 2]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 10]) + " HDD"
        preemptible: select_first([preepmtible_attempts, 2])
    }
    output {
        File lego_plots = "${output_prefix}.pdf"
    }
}

workflow legoPlotter {
    File input_file
    String input_file_format
    String output_prefix
    String sample_id
    Boolean is_exome

    String plotter_docker
    File? plotter_override
    File? renderer_override
    File? callable_regions
    File? ref_fasta

    call generateLegoPlot {
        input: plotter_docker = plotter_docker,
               plotter_override = plotter_override,
               renderer_override = renderer_override,
               input_file = input_file,
               input_file_format = input_file_format,
               output_prefix = output_prefix,
               sample_id = sample_id,
               is_exome = is_exome,
               callable_regions = callable_regions,
               ref_fasta = ref_fasta
    }

    output {
        File lego_plots = generateLegoPlot.lego_plots
    }
}
