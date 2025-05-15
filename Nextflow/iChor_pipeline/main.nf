#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

//-----------------------------------------------------
//              parameters (with sane defaults)
//-----------------------------------------------------
params.sample_id                = 'SAMPLE'          // overwrite on the CLI
params.bam_file                 = null              //  "
params.bam_index                = null              //  "
params.bin_size                 = 1000000           // 1 Mb bins

/* ---------- ichorCNA-specific ------------ */
params.gcWig                    = 'gs://gptag/ichorCNA_resources/gc_hg19_1000kb.wig'
params.mapWig                   = 'gs://gptag/ichorCNA_resources/map_hg19_1000kb.wig'
params.centromere               = 'gs://gptag/ichorCNA_resources/GRCh37.p13_centromere_UCSC-gapTable.txt'
params.normalPanel              = 'gs://gptag/ichorCNA_resources/Stover_Lennon_UH2_20HD_PON_1X_ULP_median.rds'
params.genomeBuild              = 'hg19'
params.genomeStyle              = 'NCBI'
params.chrs                     = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22'
params.chrTrain                 = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,21,22'
params.chrNormalize             = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22'
params.ploidy                   = '2,3,4'
params.normal                   = 'c(0.35,0.45,0.55,0.65,0.75,0.85,0.95)'
params.maxCN                    = 7
params.includeHOMD              = false
params.estimateNormal           = true
params.estimatePloidy           = true
params.estimateClonality        = true
params.scStates                 = 'c(1,3)'
params.txnE                     = 0.99999
params.txnStrength              = 100000
params.minSegmentBins           = 25
params.minMapScore              = 0.9
params.fracReadsChrYMale        = 0.002
params.maxFracCNASubclone       = 0.7
params.maxFracGenomeSubclone    = 0.5
params.altFracThreshold         = 0.7
params.lambdaScaleHyperParam    = 3
// any other params from the JSON may be appended here
//-----------------------------------------------------

//-----------------------------------------------------
//                 process: read_counter
//-----------------------------------------------------
process read_counter {

    /* publishDir intentionally removed (caused trouble in ICA) */

    tag { sample_id }
    container '079623148045.dkr.ecr.us-east-1.amazonaws.com/cp-prod/506be209-106a-45b1-9f9e-a203cad509d9:latest'
    cpus 1
    memory '4 GB'
    errorStrategy 'retry'
    maxRetries 1

    input:
    path bam_file
    path bam_index
    val  sample_id
    val  chrs
    val  bin_size

    output:
    path "${sample_id}.bin${bin_size}.wig", emit: wig_file

    script:
    """
    ln -vs ${bam_index} ${sample_id}.bam.bai
    ln -vs ${bam_file}  ${sample_id}.bam
    /HMMcopy/bin/readCounter ${sample_id}.bam \
        -c '${chrs}' -w ${bin_size} -q ${params.min_qual ?: 0} \
        > ${sample_id}.bin${bin_size}.wig
    """
}

//-----------------------------------------------------
//                 process: ichor_cna
//-----------------------------------------------------
process ichor_cna {
    tag { sample_id }
    container '079623148045.dkr.ecr.us-east-1.amazonaws.com/cp-prod/b9c33152-acf1-4294-9f8e-ba2ba22da04d:latest'
    cpus 4
    memory '16 GB'
    errorStrategy 'retry'
    maxRetries 1

    input:
    path wig_file                                    // from read_counter
    val  sample_id
    val  input_gcWig
    val  input_mapWig
    val  centromere
    val  normalPanel
    val  genomeBuild
    val  genomeStyle
    val  chrTrain
    val  chrNormalize
    val  ploidy
    val  normal
    val  maxCN
    val  includeHOMD
    val  estimateNormal
    val  estimatePloidy
    val  estimateClonality
    val  scStates
    val  txnE
    val  txnStrength
    val  minSegmentBins
    val  minMapScore
    val  fracReadsChrYMale
    val  maxFracCNASubclone
    val  maxFracGenomeSubclone
    val  altFracThreshold
    val  lambdaScaleHyperParam

    output:
    path "${sample_id}_segments.txt", emit: segments_txt

    script:
    """
        Rscript /runIchorCNA.R --id "${sample_id}" \
            --outDir ./ --libdir /ichorCNA \
            --WIG "${wig_file}" \
            --gcWig  "${input_gcWig}" \
            --mapWig "${input_mapWig}" \
            --normalPanel "${normalPanel}" \
            --ploidy "${ploidy}" \
            --normal "${normal}" \
            --maxCN ${maxCN} \
            --includeHOMD ${includeHOMD} \
            --chrs "${chrs}" \
            --chrTrain "${chrTrain}" \
            --chrNormalize "${chrNormalize}" \
            --genomeStyle "${genomeStyle}" \
            --genomeBuild "${genomeBuild}" \
            --estimateNormal ${estimateNormal} \
            --estimatePloidy ${estimatePloidy}  \
            --estimateScPrevalence ${estimateClonality} \
            --scStates "${scStates}" \
            --centromere ${centromere} \
            --exons.bed ${exons} \
            --txnE ${txnE} \
            --txnStrength ${txnStrength} \
            --minSegmentBins ${minSegmentBins} \
            --minMapScore ${minMapScore} \
            --lambdaScaleHyperParam ${lambdaScaleHyperParam} \
            --fracReadsInChrYForMale ${fracReadsChrYMale} \
            --maxFracGenomeSubclone ${maxFracGenomeSubclone} \
            --altFracThreshold ${altFracThreshold} \
            --maxFracCNASubclone ${maxFracCNASubclone} \
            --rmCentromereFlankLength ${rmCentromereFlankLength} \
            --plotFileType ${plotFileType} \
            --plotYLim "${plotYlim}"
    """
}

//-----------------------------------------------------
//                    workflow
//-----------------------------------------------------
workflow {

    /* ---------- channels for read_counter ---------- */
    Channel.fromPath(params.bam_file)  .set { bam_file_ch  }
    Channel.fromPath(params.bam_index) .set { bam_index_ch }
    Channel.value(params.sample_id)    .set { sample_id_ch }
    Channel.value(params.chrs)         .set { chrs_ch }
    Channel.value(params.bin_size)     .set { bin_size_ch }

    /* run read_counter */
    read_counter_out = read_counter(
        bam_file_ch,
        bam_index_ch,
        sample_id_ch,
        chrs_ch,
        bin_size_ch
    )

    /* ---------- channels for ichor_cna ---------- */
    Channel.value(params.gcWig)                 .set { gcWig_ch }
    Channel.value(params.mapWig)                .set { mapWig_ch }
    Channel.value(params.centromere)            .set { centromere_ch }
    Channel.value(params.normalPanel)           .set { normalPanel_ch }
    Channel.value(params.genomeBuild)           .set { genomeBuild_ch }
    Channel.value(params.genomeStyle)           .set { genomeStyle_ch }
    Channel.value(params.chrTrain)              .set { chrTrain_ch }
    Channel.value(params.chrNormalize)          .set { chrNormalize_ch }
    Channel.value(params.ploidy)                .set { ploidy_ch }
    Channel.value(params.normal)                .set { normal_ch }
    Channel.value(params.maxCN)                 .set { maxCN_ch }
    Channel.value(params.includeHOMD)           .set { includeHOMD_ch }
    Channel.value(params.estimateNormal)        .set { estimateNormal_ch }
    Channel.value(params.estimatePloidy)        .set { estimatePloidy_ch }
    Channel.value(params.estimateClonality)     .set { estimateClonality_ch }
    Channel.value(params.scStates)              .set { scStates_ch }
    Channel.value(params.txnE)                  .set { txnE_ch }
    Channel.value(params.txnStrength)           .set { txnStrength_ch }
    Channel.value(params.minSegmentBins)        .set { minSegmentBins_ch }
    Channel.value(params.minMapScore)           .set { minMapScore_ch }
    Channel.value(params.fracReadsChrYMale)     .set { fracReadsChrYMale_ch }
    Channel.value(params.maxFracCNASubclone)    .set { maxFracCNASubclone_ch }
    Channel.value(params.maxFracGenomeSubclone) .set { maxFracGenomeSubclone_ch }
    Channel.value(params.altFracThreshold)      .set { altFracThreshold_ch }
    Channel.value(params.lambdaScaleHyperParam) .set { lambdaScaleHyperParam_ch }

    /* run ichor_cna */
    ichor_cna(
        read_counter_out.wig_file,
        sample_id_ch,
        gcWig_ch,
        mapWig_ch,
        centromere_ch,
        normalPanel_ch,
        genomeBuild_ch,
        genomeStyle_ch,
        chrTrain_ch,
        chrNormalize_ch,
        ploidy_ch,
        normal_ch,
        maxCN_ch,
        includeHOMD_ch,
        estimateNormal_ch,
        estimatePloidy_ch,
        estimateClonality_ch,
        scStates_ch,
        txnE_ch,
        txnStrength_ch,
        minSegmentBins_ch,
        minMapScore_ch,
        fracReadsChrYMale_ch,
        maxFracCNASubclone_ch,
        maxFracGenomeSubclone_ch,
        altFracThreshold_ch,
        lambdaScaleHyperParam_ch
    )
}

