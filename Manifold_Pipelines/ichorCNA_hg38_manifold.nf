#!/usr/bin/env nextflow
/*
 * ichorCNA (hg38) — Nextflow DSL2 port of iChorCNA_hg38.wdl
 *
 * Pipeline:
 *   READ_COUNTER  ->  ICHORCNA  ->  EXTRACT_ICHOR_PARAMS
 *                             \->  BUNDLE_PER_CHROMOSOME_PLOTS
 *
 * Run:
 *   nextflow run ichorCNA_hg38.nf -params-file ichorCNA_hg38.params.json
 */

nextflow.enable.dsl = 2

/* -------------------------------------------------------------------------- *
 *  Processes
 * -------------------------------------------------------------------------- */

process READ_COUNTER {
    tag    { sample_id }
    label  'read_counter'
    container params.read_counter_docker

    input:
    tuple val(sample_id), path(bam_file), path(bam_index)

    output:
    tuple val(sample_id), path("${sample_id}.bin${params.bin_size}.wig"), emit: wig

    script:
    """
    set -euxo pipefail

    ln -vs "${bam_index}" "${sample_id}.bam.bai"
    ln -vs "${bam_file}"  "${sample_id}.bam"

    /HMMcopy/bin/readCounter ${sample_id}.bam \\
        -c ${params.chr_counter} \\
        -w ${params.bin_size} \\
        -q ${params.read_counter_min_qual} \\
        > "${sample_id}.bin${params.bin_size}.wig"
    """
}

process ICHORCNA {
    tag    { sample_id }
    label  'ichorcna'
    container params.ichorcna_docker
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(wig_file)
    path  gcWig
    path  mapWig
    path  centromere
    path  normalPanel          // stageAs 'NO_FILE' when absent (see workflow)
    path  exons                // stageAs 'NO_FILE' when absent (see workflow)

    output:
    tuple val(sample_id), path("${sample_id}.params.txt"),                        emit: params
    tuple val(sample_id), path("${sample_id}.correctedDepth.txt"),                emit: corrDepth
    tuple val(sample_id), path("${sample_id}.cna.seg"),                           emit: cna
    tuple val(sample_id), path("${sample_id}.seg.txt"),                           emit: segTxt
    tuple val(sample_id), path("${sample_id}.seg"),                               emit: seg
    tuple val(sample_id), path("${sample_id}.RData"),                             emit: rdata
    tuple val(sample_id), path("${sample_id}/${sample_id}_genomeWide_all_sols.pdf"), emit: allGenomeWidePlots
    tuple val(sample_id), path("${sample_id}/${sample_id}_bias.pdf"),             emit: bias
    tuple val(sample_id), path("${sample_id}/${sample_id}_tpdf.pdf"),             emit: tpdf
    tuple val(sample_id), path("${sample_id}/${sample_id}_correct.pdf"),          emit: correct
    tuple val(sample_id), path("${sample_id}/${sample_id}_CNA*"),                 emit: perChromosomePlots
    tuple val(sample_id), path("${sample_id}.optimalSolution.zip"),               emit: optimalSolution
    tuple val(sample_id), path("${sample_id}.outSolutions.zip"),                  emit: outSolutions

    script:
    // Optional-file handling: an empty stub named NO_FILE means "not provided"
    def normalPanelArg = normalPanel.name != 'NO_FILE' ? normalPanel : 'None'
    def exonsArg       = exons.name       != 'NO_FILE' ? "--exons.bed '${exons}'" : ''
    def coverageArg    = params.mean_depth != null ? "--coverage ${params.mean_depth}" : ''
    """
    set -euxo pipefail

    Rscript /runIchorCNA.R --id "${sample_id}" \\
        --outDir ./ --libdir /ichorCNA \\
        --WIG "${wig_file}" \\
        --gcWig  "${gcWig}" \\
        --mapWig "${mapWig}" \\
        --normalPanel "${normalPanelArg}" \\
        --ploidy "${params.ploidy}" \\
        --normal "${params.normal}" \\
        ${coverageArg} \\
        --maxCN ${params.maxCN} \\
        --includeHOMD ${params.includeHOMD} \\
        --chrs "${params.chrs}" \\
        --chrTrain "${params.chrTrain}" \\
        --chrNormalize "${params.chrNormalize}" \\
        --genomeStyle "${params.genome_style}" \\
        --genomeBuild "${params.genomeBuild}" \\
        --estimateNormal ${params.estimateNormal} \\
        --estimatePloidy ${params.estimatePloidy} \\
        --estimateScPrevalence ${params.estimateClonality} \\
        --scStates "${params.scStates}" \\
        --centromere ${centromere} \\
        ${exonsArg} \\
        --txnE ${params.txnE} \\
        --txnStrength ${params.txnStrength} \\
        --minSegmentBins ${params.minSegmentBins} \\
        --minMapScore ${params.minMapScore} \\
        --lambdaScaleHyperParam ${params.lambdaScaleHyperParam} \\
        --fracReadsInChrYForMale ${params.fracReadsChrYMale} \\
        --maxFracGenomeSubclone ${params.maxFracGenomeSubclone} \\
        --altFracThreshold ${params.altFracThreshold} \\
        --maxFracCNASubclone ${params.maxFracCNASubclone} \\
        --rmCentromereFlankLength ${params.rmCentromereFlankLength} \\
        --plotFileType ${params.plotFileType} \\
        --plotYLim "${params.plotYlim}"

    # Zip optimal solutions
    mkdir "${sample_id}.optimalSolution"
    cp "${sample_id}/${sample_id}_genomeWide.pdf" \\
        "${sample_id}.cna.seg" \\
        "${sample_id}.seg.txt" \\
        "${sample_id}.seg" \\
        "${sample_id}.optimalSolution/"
    zip -r "${sample_id}.optimalSolution.zip" "${sample_id}.optimalSolution"

    # Generate list of out solutions
    Rscript /gatherOutSolutions.R --id "${sample_id}" --rdata "${sample_id}.RData"
    """
}

process EXTRACT_ICHOR_PARAMS {
    tag    { sample_id }
    label  'extract_params'
    container params.python_docker
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(params_file)

    output:
    tuple val(sample_id),
          path("tumor_fraction"),
          path("coverage"),
          path("ploidy"),
          path("subclone_fraction"),
          path("fraction_genome_subclonal"),
          path("fraction_cna_subclonal"),
          path("gc-map_correction_mad"),
          path("top_solution_log_likelihood"), emit: values

    script:
    """
    set -euxo pipefail

    cut -f1,2 "${params_file}" | grep -v ^\$ | tr "\\t" " " > params_table.txt
python3<<CODE
with open("${params_file}", "r") as solutions_in:
    solutions = solutions_in.readlines()
sols = [x.rstrip("\\n").split("\\t") for x in solutions if x.startswith("n0.")]
log_lik = 0
for sol in sols:
    if (
        int(float(sol[6])) > log_lik
        and float(sol[4]) < ${params.maxFracGenomeSubclone}
        and float(sol[5]) < ${params.maxFracCNASubclone}
    ):
        log_lik = int(float(sol[6]))
params = open("params_table.txt", "r").readlines()
params = [x.rstrip("\\n").split(": ") for x in params if ":" in x]
param_dict = {a: b for a, b in params}
with open("tumor_fraction", "w") as p:
    p.write(param_dict["Tumor Fraction"])
with open("coverage", "w") as p:
    p.write(param_dict["Coverage"])
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
    """
}

process BUNDLE_PER_CHROMOSOME_PLOTS {
    tag    { sample_id }
    label  'bundle_plots'
    container params.tagtools_docker
    publishDir "${params.outdir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(chrom_plots)

    output:
    tuple val(sample_id), path("${sample_id}_OptimalSolutionPerChrom.pdf"), emit: output_plot

    script:
    """
    set -euxo pipefail

    # load list of plots into an array, which will correctly handle files with spaces in them
    printf '%s\\n' ${chrom_plots} | sort -V > chrom_plots.txt
    readarray -t CHROM_PLOTS < chrom_plots.txt

    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite \\
        -sOutputFile="${sample_id}_OptimalSolutionPerChrom.pdf" \\
        "\${CHROM_PLOTS[@]}"
    """
}

/* -------------------------------------------------------------------------- *
 *  Workflow
 * -------------------------------------------------------------------------- */

workflow {

    // --- inputs ---
    ch_bam = Channel.of( tuple(
        params.sample_id,
        file(params.bam_file),
        file(params.bam_index)
    ) )

    // Optional files: fall back to an empty stub named NO_FILE
    def no_file    = file("${projectDir}/assets/NO_FILE")
    ch_gcWig       = file(params.gcWig)
    ch_mapWig      = file(params.mapWig)
    ch_centromere  = file(params.centromere)
    ch_normalPanel = params.normalPanel ? file(params.normalPanel) : no_file
    ch_exons       = params.exons       ? file(params.exons)       : no_file

    // --- run ---
    READ_COUNTER( ch_bam )

    ICHORCNA(
        READ_COUNTER.out.wig,
        ch_gcWig,
        ch_mapWig,
        ch_centromere,
        ch_normalPanel,
        ch_exons
    )

    EXTRACT_ICHOR_PARAMS( ICHORCNA.out.params )

    BUNDLE_PER_CHROMOSOME_PLOTS( ICHORCNA.out.perChromosomePlots )

    // --- QC gate (mirrors WDL `qc_passed`) ---
    EXTRACT_ICHOR_PARAMS.out.values
        .map { sid, tf, cov, plo, scf, fgs, fcs, mad, ll ->
            def coverage = cov.text.trim() as Float
            def madScore = mad.text.trim() as Float
            def passed   = (madScore <= params.max_passing_mad_score) &&
                           (coverage >= params.min_passing_coverage)
            tuple(sid, passed, coverage, madScore)
        }
        .subscribe { sid, passed, coverage, madScore ->
            log.info "[${sid}] qc_passed=${passed} (coverage=${coverage}, gc_map_mad=${madScore})"
        }
}
