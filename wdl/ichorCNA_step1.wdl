#upload commands
#docker run -it --rm -v "$HOME"/.config:/.config -v "$PWD":/working broadinstitute/firecloud-cli /bin/bash
#firecloud -m push -s chrislo -n ichorcna -t Workflow -y "ichorcna" ichorcna.wdl 

#todo:
#put back NORMWIG, lambda etc.
#put back "X"

workflow ichorCNA {
   # Disk size parameter
   Int disk_size

   # bamToWigTask
   File bam_file
   File bam_index_file
   String sample_id

   # readsCorrectionTask
   Boolean normalizeMaleX
   String chrs
   String chrNormalize
   File centromere
   File normalPanel
   #File targetedSequence
   #File NORMWIG
   Int rmCentromereFlankLength

   # HMMTask
   String ploidyParams
   Boolean diploidChrX
   String normal
   Int maxCN
   Boolean estimateNormal
   Boolean estimatePloidy
   Boolean includeHOMD
   String chrTrain
   Float txnE
   Int txnStrength
   String scStates
   Boolean estimateScPrevalence
   Float maxFracGenomeSubclone
   Float maxFracCNASubclone
   Int minSegmentBins
   Float altFracThreshold
   Int lambdaScaleHyperParam
   Float mean_depth

   call bamToWigTask {
	   input: disk_size = disk_size,
              bam_file = bam_file,
              bam_index_file = bam_index_file,
              sample_id = sample_id
   }

   call readsCorrectionTask {
      input: disk_size = disk_size,
             wig_file = bamToWigTask.wig_file,
             sample_id = sample_id,
             normalizeMaleX = normalizeMaleX,
             chrs = chrs,
             chrNormalize = chrNormalize,
             centromere = centromere,
             normalPanel = normalPanel,
             rmCentromereFlankLength = rmCentromereFlankLength
   }

   call HMMTask {
      input: disk_size = disk_size,
             tumorCorrectedDepth_Rdata = readsCorrectionTask.tumorCorrectedDepthRdata,
             sample_id = sample_id,
             ploidyParams = ploidyParams,
             diploidChrX = diploidChrX,
             normal = normal,
             maxCN = maxCN,
             estimateNormal = estimateNormal,
             estimatePloidy = estimatePloidy,
             includeHOMD = includeHOMD,
             chrs = chrs,
             chrTrain = chrTrain,
             txnE = txnE,
             txnStrength = txnStrength,
             scStates = scStates,
             estimateScPrevalence = estimateScPrevalence,
             maxFracGenomeSubclone = maxFracGenomeSubclone,
             maxFracCNASubclone = maxFracCNASubclone,
             minSegmentBins = minSegmentBins,
             altFracThreshold = altFracThreshold,
             lambdaScaleHyperParam = lambdaScaleHyperParam,
             mean_depth = mean_depth
   }

   call BundlePerChromosomePlots {
       input: chrom_plots = HMMTask.perChromosomePlots,
              sample_id = sample_id
   }

   output {
      File allGenomeWidePlots = HMMTask.allGenomeWidePlots

      Float tumor_fraction = HMMTask.tumor_fraction
      Float ploidy = HMMTask.ploidy
      String subclone_fraction = HMMTask.subclone_fraction
      String fraction_genome_subclonal = HMMTask.fraction_genome_subclonal
      String fraction_cna_subclonal = HMMTask.fraction_cna_subclonal
      String gc_map_correction_mad = HMMTask.gc_map_correction_mad

      File bias = HMMTask.bias
      File tpdf = HMMTask.tpdf
      File correct = HMMTask.correct
      File params = HMMTask.params
      
      File optimalSolution = HMMTask.optimalSolution
      File outSolutions = HMMTask.outSolutions
      File perChromosomePlots = BundlePerChromosomePlots.output_plot
   }
}


task bamToWigTask {
   Int disk_size
   File bam_file
   File bam_index_file
   String sample_id

   command <<<
      ln -vs ${bam_index_file} ${sample_id}.bam.bai
      ln -vs ${bam_file} ${sample_id}.bam

      /HMMcopy/bin/readCounter --window 1000000 \
                               --quality 20 \
                               --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" \
                               ${sample_id}.bam > ${sample_id}.wig
   >>>
   output {
      File wig_file = "${sample_id}.wig"
   }
   runtime {
      docker: "us.gcr.io/tag-team-160914/bloodbiopsy-hmmcopy:0.0.1"
      disks: "local-disk ${disk_size} HDD"
      preemptible: 3
   }
}


task readsCorrectionTask {
   Int disk_size
   File wig_file
   String sample_id
   Boolean normalizeMaleX
   String chrs
   String chrNormalize
   File centromere
   File normalPanel
   #File targetedSequence
   #File NORMWIG
   Int rmCentromereFlankLength

   command <<<
      Rscript /readsCorrection.R --libdir /ichorCNA/R/ \
                                 --datadir /ichorCNA/inst/extdata/ \
                                 --id ${sample_id} \
                                 --WIG ${wig_file} \
                                 --normalizeMaleX ${normalizeMaleX} \
                                 --chrs "${chrs}" \
                                 --chrNormalize "${chrNormalize}" \
                                 --centromere ${centromere} \
                                 --normalPanel ${normalPanel} \
                                 --rmCentromereFlankLength ${rmCentromereFlankLength} \
                                 --outDir . 
   >>>
   output {
      File correctedDepth = "${sample_id}/${sample_id}.correctedDepth.txt"
      File tumorCorrectedDepthRdata = "${sample_id}/${sample_id}.tumor_copy_correctedDepth.RData"
   }
   runtime {
      docker: "us.gcr.io/tag-team-160914/bloodbiopsy-ichorcna:0.0.1"
      disks: "local-disk ${disk_size} HDD"
      preemptible: 3
   }
}


task HMMTask {
   Int disk_size
   File tumorCorrectedDepth_Rdata
   String sample_id
   String ploidyParams
   Boolean diploidChrX
   String normal
   Int maxCN
   Boolean estimateNormal
   Boolean estimatePloidy
   Boolean includeHOMD
   String chrs
   String chrTrain
   Float txnE
   Int txnStrength
   String scStates
   Boolean estimateScPrevalence
   Float maxFracGenomeSubclone
   Float maxFracCNASubclone
   Int minSegmentBins
   Float altFracThreshold
   Int lambdaScaleHyperParam
   Float mean_depth

   command <<<
      Rscript /HMM.R --libdir /ichorCNA/R/ \
                     --id ${sample_id} \
                     --tumourCopyCorrected ${tumorCorrectedDepth_Rdata} \
                     --ploidy "${ploidyParams}" \
                     --diploidChrX ${diploidChrX} \
                     --normal "${normal}" \
                     --maxCN ${maxCN} \
                     --estimateNormal ${estimateNormal} \
                     --estimatePloidy ${estimatePloidy} \
                     --includeHOMD ${includeHOMD} \
                     --chrs "${chrs}" \
                     --chrTrain "${chrTrain}" \
                     --txnE ${txnE} \
                     --txnStrength ${txnStrength} \
                     --scStates "${scStates}" \
                     --estimateScPrevalence ${estimateScPrevalence} \
                     --maxFracGenomeSubclone ${maxFracGenomeSubclone} \
                     --maxFracCNASubclone ${maxFracCNASubclone} \
                     --minSegmentBins ${minSegmentBins} \
                     --altFracThreshold ${altFracThreshold} \
                     --lambdaScaleHyperParam ${lambdaScaleHyperParam} \
                     --coverage ${mean_depth} --outDir .

      # Extract sample stats: though a bit ugly, grep'ing all the values
      grep "^Tumor Fraction" ${sample_id}/${sample_id}.params.txt | cut -f2 > tumor_fraction
      grep "^Ploidy" ${sample_id}/${sample_id}.params.txt | cut -f2 > ploidy
      grep "^Subclone Fraction" ${sample_id}/${sample_id}.params.txt | cut -f2 > subclone_fraction
      grep "^Fraction Genome Subclonal" ${sample_id}/${sample_id}.params.txt | cut -f2 > fraction_genome_subclonal
      grep "^Fraction CNA Subclonal" ${sample_id}/${sample_id}.params.txt | cut -f2 > fraction_cna_subclonal
      grep "^GC-Map correction MAD" ${sample_id}/${sample_id}.params.txt | cut -f2 > gc-map_correction_mad

      # Zip optimal solutions
      mkdir ${sample_id}.optimalSolution
      mv ${sample_id}/optimalSolution/* ${sample_id}.optimalSolution
      zip -r ${sample_id}.optimalSolution.zip ${sample_id}.optimalSolution
      
      # Zip the other solutions
      mkdir ${sample_id}.outSolutions
      mv ${sample_id}/n*/* ${sample_id}.outSolutions
      zip -r ${sample_id}.outSolutions.zip ${sample_id}.outSolutions
   >>>
   output {
      File allGenomeWidePlots = "${sample_id}/${sample_id}_genomeWide_all_sols.pdf"

      Float tumor_fraction = read_float("tumor_fraction")
      Float ploidy = read_float("ploidy")
      String subclone_fraction = read_string("subclone_fraction")
      String fraction_genome_subclonal = read_string("fraction_genome_subclonal")
      String fraction_cna_subclonal = read_string("fraction_cna_subclonal")
      String gc_map_correction_mad = read_string("gc-map_correction_mad")

      File bias = "${sample_id}/${sample_id}_bias.pdf"
      File tpdf = "${sample_id}/${sample_id}_tpdf.pdf"
      File correct = "${sample_id}/${sample_id}_correct.pdf"
      File params = "${sample_id}/${sample_id}.params.txt"
      
      File optimalSolution = "${sample_id}.optimalSolution.zip"
      File outSolutions = "${sample_id}.outSolutions.zip"
      Array[File] perChromosomePlots = glob("${sample_id}/${sample_id}_CNA*")
   }
   runtime {
      docker: "us.gcr.io/tag-team-160914/bloodbiopsy-ichorcna:0.0.1"
      disks: "local-disk ${disk_size} HDD"
      preemptible: 3
   }
}

task BundlePerChromosomePlots {
   Array[File] chrom_plots
   String sample_id

   command <<<
      set -e
      CHROM_PLOTS=`ls ${sep=" " chrom_plots} | sort -V`

      gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite \
         -sOutputFile="${sample_id}_OptimalSolutionPerChrom.pdf" \
         $CHROM_PLOTS
   >>>
   output {
      File output_plot = "${sample_id}_OptimalSolutionPerChrom.pdf"
   }
   runtime {
      docker: "us.gcr.io/tag-team-160914/tag-tools:0.0.4"
      memory: "1 GB"
      disks: "local-disk 2 HDD"
      preemptible: 3
   }
}
