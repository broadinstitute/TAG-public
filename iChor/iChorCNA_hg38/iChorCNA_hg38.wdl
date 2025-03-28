version 1.0

workflow ichorCNA {
    input {
        File bam_file
        File bam_index
        String sample_id
        String genomeBuild
        File? normalPanel
        Int bin_size_kb
        File centromere
        File gcWig
        File mapWig
        String genome_style
        String chr_counter # All chromosomes for counting reads, as a list
        String chrs # Autosomal and female chromosome for running ichor, as an R command
        String chrTrain # Autosomal chromosomes to estimate ichor params, as an R command
        Float maxFracCNASubclone
        Float maxFracGenomeSubclone
    }

    Int bin_size = bin_size_kb * 1000

    call read_counter {
        input:
            bam_file = bam_file,
            bam_index = bam_index,
            sample_id = sample_id,
            bin_size = bin_size,
            chrs = chr_counter,
    }

    call ichorCNATask {
        input:
            wig_file = read_counter.wig_file,
            genomeBuild = genomeBuild,
            sample_id = sample_id,
            normalPanel = normalPanel,
            chrs = chrs,
            chrTrain = chrTrain,
            centromere = centromere,
            genomeStyle = genome_style,
            input_gcWig = gcWig,
            input_mapWig = mapWig,
            maxFracCNASubclone = maxFracCNASubclone,
            maxFracGenomeSubclone = maxFracGenomeSubclone,
    }
    call extractIchorParams {
        input:
            params = ichorCNATask.params,
            maxFracCNASubclone = maxFracCNASubclone,
            maxFracGenomeSubclone = maxFracGenomeSubclone,
    }

    call bundlePerChromosomePlots {
        input:
            chrom_plots = ichorCNATask.perChromosomePlots,
            sample_id = sample_id,
    }

    output {
        File wig_file = read_counter.wig_file

        File allGenomeWidePlots = ichorCNATask.allGenomeWidePlots
        File corrDepth = ichorCNATask.corrDepth
        File cna = ichorCNATask.cna
        File segTxt = ichorCNATask.segTxt
        File seg = ichorCNATask.seg
        File rdata = ichorCNATask.rdata

        Float tumor_fraction = extractIchorParams.tumor_fraction
        Float ploidy = extractIchorParams.ploidy
        String subclone_fraction = extractIchorParams.subclone_fraction
        String fraction_genome_subclonal = extractIchorParams.fraction_genome_subclonal
        String fraction_cna_subclonal = extractIchorParams.fraction_cna_subclonal
        String gc_map_correction_mad = extractIchorParams.gc_map_correction_mad
        Int top_solution_log_likelihood = extractIchorParams.top_solution_log_likelihood

        File bias = ichorCNATask.bias
        File tpdf = ichorCNATask.tpdf
        File correct = ichorCNATask.correct
        File params = ichorCNATask.params

        File optimalSolution = ichorCNATask.optimalSolution
        File outSolutions = ichorCNATask.outSolutions
        File perChromosomePlots = bundlePerChromosomePlots.output_plot
    }
}

task read_counter {
    input {
        File bam_file
        File bam_index
        String chrs
        String sample_id
        Int bin_size

        Int min_qual = 20
        Int memGB = 4
        Int diskGB = 50
        Int preemptible = 2

        String docker_override = "us.gcr.io/tag-public/bloodbiopsy-hmmcopy:0.0.1"
    }

    String wig_file_name = "~{sample_id}.bin~{bin_size}.wig"

    command <<<
        set -euxo pipefail

        ln -vs "~{bam_index}" "~{sample_id}.bam.bai"
        ln -vs "~{bam_file}" "~{sample_id}.bam"

        /HMMcopy/bin/readCounter ~{sample_id}.bam -c ~{chrs} -w ~{bin_size} -q ~{min_qual} \
            > "~{wig_file_name}"
    >>>

    runtime {
		docker: select_first([docker_override, "us.gcr.io/tag-public/bloodbiopsy-hmmcopy:0.0.1"])
        disks: "local-disk ~{diskGB} HDD"
        memory: "~{memGB}GB"
        preemptible: "~{preemptible}"
        maxRetries: 1
    }

    output {
        File wig_file = wig_file_name
    }
}

task ichorCNATask {
    input {
        File wig_file
        File? normalPanel
        String genomeBuild
        Float? mean_depth
        String sample_id
        String ploidy
        String normal
        Int maxCN
        Boolean includeHOMD
        String chrs
        String chrTrain
        String chrNormalize
        String genomeStyle
        Boolean estimateNormal
        Boolean estimatePloidy
        Boolean estimateClonality
        String scStates
        File centromere
        File? exons
        Float txnE
        Int txnStrength
        Int minSegmentBins
        Float minMapScore
        Float fracReadsChrYMale
        Float maxFracCNASubclone
        Float maxFracGenomeSubclone
        Float altFracThreshold
        Int lambdaScaleHyperParam
        Int rmCentromereFlankLength
        String plotFileType = "pdf"
        String plotYlim
        File input_gcWig
        File input_mapWig

        Int memGB = 4
        Int diskGB = 50
        Int preemptible = 2
        Int maxRetries = 1

        String docker_override = "us.gcr.io/tag-public/bloodbiopsy-ichorcna:0.2.1"
    }

    command <<<
        set -euxo pipefail

        Rscript /runIchorCNA.R --id "~{sample_id}" \
            --outDir ./ --libdir /ichorCNA \
            --WIG "~{wig_file}" \
            --gcWig  "~{input_gcWig}" \
            --mapWig "~{input_mapWig}" \
            --normalPanel "~{default="None" normalPanel}" \
            --ploidy "~{ploidy}" \
            --normal "~{normal}" \
            ~{"--coverage " + mean_depth} \
            --maxCN ~{maxCN} \
            --includeHOMD ~{includeHOMD} \
            --chrs "~{chrs}" \
            --chrTrain "~{chrTrain}" \
            --chrNormalize "~{chrNormalize}" \
            --genomeStyle "~{genomeStyle}" \
            --genomeBuild "~{genomeBuild}" \
            --estimateNormal ~{estimateNormal} \
            --estimatePloidy ~{estimatePloidy}  \
            --estimateScPrevalence ~{estimateClonality} \
            --scStates "~{scStates}" \
            --centromere ~{centromere} \
            ~{"--exons.bed '" + exons + "'"} \
            --txnE ~{txnE} \
            --txnStrength ~{txnStrength} \
            --minSegmentBins ~{minSegmentBins} \
            --minMapScore ~{minMapScore} \
            --lambdaScaleHyperParam ~{lambdaScaleHyperParam} \
            --fracReadsInChrYForMale ~{fracReadsChrYMale} \
            --maxFracGenomeSubclone ~{maxFracGenomeSubclone} \
            --altFracThreshold ~{altFracThreshold} \
            --maxFracCNASubclone ~{maxFracCNASubclone} \
            --rmCentromereFlankLength ~{rmCentromereFlankLength} \
            --plotFileType ~{plotFileType} \
            --plotYLim "~{plotYlim}"

        # Zip optimal solutions
        mkdir "~{sample_id}.optimalSolution"
        cp "~{sample_id}/~{sample_id}_genomeWide.pdf" \
            "~{sample_id}.cna.seg" \
            "~{sample_id}.seg.txt" \
            "~{sample_id}.seg" \
            "~{sample_id}.optimalSolution/"
        zip -r "~{sample_id}.optimalSolution.zip" "~{sample_id}.optimalSolution"

        # Generate list of out solutions
        Rscript /gatherOutSolutions.R --id "~{sample_id}" --rdata "~{sample_id}.RData"
    >>>

    runtime {
		docker: select_first([docker_override, "us.gcr.io/tag-public/bloodbiopsy-ichorcna:0.2.1"])
        disks: "local-disk ~{diskGB} HDD"
        memory: "~{memGB} GB"
        preemptible: preemptible
        maxRetries: maxRetries
    }

    output {
        File corrDepth = "~{sample_id}.correctedDepth.txt"
        File params = "~{sample_id}.params.txt"
        File cna = "~{sample_id}.cna.seg"
        File segTxt = "~{sample_id}.seg.txt"
        File seg = "~{sample_id}.seg"
        File rdata = "~{sample_id}.RData"

        File allGenomeWidePlots = "~{sample_id}/~{sample_id}_genomeWide_all_sols.pdf"
        File bias = "~{sample_id}/~{sample_id}_bias.pdf"
        File tpdf = "~{sample_id}/~{sample_id}_tpdf.pdf"
        File correct = "~{sample_id}/~{sample_id}_correct.pdf"
        Array[File] perChromosomePlots = glob("~{sample_id}/~{sample_id}_CNA*")

        File optimalSolution = "~{sample_id}.optimalSolution.zip"
        File outSolutions = "~{sample_id}.outSolutions.zip"
    }

}


task extractIchorParams {
    input {
        File params
        Float maxFracCNASubclone
        Float maxFracGenomeSubclone
    }

    command <<<
        set -euxo pipefail

        cut -f1,2 "~{params}" | grep -v ^$ | tr "\t" " " > params_table.txt
python<<CODE
with open("~{params}", "r") as solutions_in:
    solutions = solutions_in.readlines()
sols = [x.rstrip("\n").split("\t") for x in solutions if x.startswith("n0.")]
log_lik = 0
for sol in sols:
    if (
        int(float(sol[6])) > log_lik
        and float(sol[4]) < ~{maxFracGenomeSubclone}
        and float(sol[5]) < ~{maxFracCNASubclone}
    ):
        log_lik = int(float(sol[6]))
params = open("params_table.txt", "r").readlines()
params = [x.rstrip("\n").split(": ") for x in params if ":" in x]
param_dict = {a: b for a, b in params}
with open("tumor_fraction", "w") as p:
    p.write(param_dict["Tumor Fraction"])
with open("ploidy", "w") as p:
    p.write(param_dict["Ploidy"])
with open("subclone_fraction", "w") as p:
    p.write(param_dict["Subclone Fraction"])
with open("fraction_genome_subclonal", "w") as p:
    p.write(param_dict["Fraction Genome Subclonal"])
with open("fraction_cna_subclonal", "w") as p:
    p.write(param_dict["Fraction CNA Subclonal"])
with open("gc-map_correction_mad", "w") as p:
    p.write(param_dict["GC-Map correction MAD"])
with open("top_solution_log_likelihood", "w") as p:
    p.write(str(log_lik))

CODE
    >>>

    runtime {
        docker: "python:3"
        memory: "2 GB"
        disks: "local-disk 10 HDD"
        preemptible: 1
        maxRetries: 1
    }

    output {
        Float tumor_fraction = read_float("tumor_fraction")
        Float ploidy = read_float("ploidy")
        String subclone_fraction = read_string("subclone_fraction")
        String fraction_genome_subclonal = read_string("fraction_genome_subclonal")
        String fraction_cna_subclonal = read_string("fraction_cna_subclonal")
        String gc_map_correction_mad = read_string("gc-map_correction_mad")
        Int top_solution_log_likelihood = read_int("top_solution_log_likelihood")
    }
}

task bundlePerChromosomePlots {
    input {
        Array[File] chrom_plots
        String sample_id

        Int memGB = 2
        Int diskGB = 10
        Int preemptible = 3
        Int maxRetries = 1

        String docker_override = "us.gcr.io/tag-public/tag-tools:1.0.0"
    }

    command <<<
        set -euxo pipefail

        # load list of plots into an array, which will correctly handle files with spaces in them
        readarray -t CHROM_PLOTS < <(sort -V "~{write_lines(chrom_plots)}")

        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite \
            -sOutputFile="~{sample_id}_OptimalSolutionPerChrom.pdf" \
            "${CHROM_PLOTS[@]}"
   >>>

   output {
      File output_plot = "~{sample_id}_OptimalSolutionPerChrom.pdf"
   }

   runtime {
        docker: select_first([docker_override, "us.gcr.io/tag-public/tag-tools:1.0.0"])
        memory: "~{memGB} GB"
        disks: "local-disk ~{diskGB} HDD"
        preemptible: preemptible
        maxRetries: maxRetries
   }
}
