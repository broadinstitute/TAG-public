version 1.0

workflow SigProfiler {
    input {
        Array[File] MutlistFiles
        File GenomeFasta
    }

    call mutlist_to_96_contexts {
        input:
            MutlistFiles = MutlistFiles,
            GenomeFasta = GenomeFasta
    }
    call sigprofiler_analysis {
        input: 
            MutationMetrics = mutlist_to_96_contexts.MutationMetrics

    }
    call PlotSignatures {
        input: 
            SignatureCount = sigprofiler_analysis.SignatureCount
    }

    output {
        File MutationMetrics = mutlist_to_96_contexts.MutationMetrics
        File SpectrumPlots = mutlist_to_96_contexts.SpectrumPlots
        File DecomposedSignatureProbabilities = sigprofiler_analysis.DecomposedSignatureProbabilities
        File SignatureStackedPlot = sigprofiler_analysis.SignatureStackedPlot
        File TMBPlot = sigprofiler_analysis.TMBPlot
        File SignatureCount = sigprofiler_analysis.SignatureCount
        File SignatureProportionPDF = PlotSignatures.signature_proportions_pdf
    }
}



task mutlist_to_96_contexts {
    input {
        Array[File] MutlistFiles
        File GenomeFasta
        File GenomeFastaIndex
    }

    command {
        Rscript /scripts/96_contexts_mutations.R "~{sep=' ' MutlistFiles}" ~{GenomeFasta}
    }
    
    output {
        File MutationMetrics = "trinuc_mutation_metrics.txt"
        File SpectrumPlots = "all_sample_spectrums.pdf"
    }

    runtime {
        docker: "us.gcr.io/tag-public/sigprofiler:v1"
        memory: "8 GB"
        disks: "local-disk 20 HDD"
    }
}

task sigprofiler_analysis {
    input {
        File MutationMetrics
        String OutputFolder = "SigProfiler-output"
    }

    command {
        python3 <<EOF
        import sys
        from SigProfilerMatrixGenerator import install as genInstall
        genInstall.install("GRCh38")

        from SigProfilerAssignment import Analyzer as Analyze
        Analyze.cosmic_fit(samples="~{MutationMetrics}",
                           output="~{OutputFolder}",
                           input_type="matrix",
                           genome_build="GRCh38",
                           cosmic_version=3.3)
        EOF
    }

    output {
        File DecomposedSignatureProbabilities = "~{OutputFolder}/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities.txt"
        File SignatureStackedPlot = "~{OutputFolder}/Assignment_Solution/Activities/Assignment_Solution_Activity_Plots.pdf"
        File TMBPlot = "~{OutputFolder}/Assignment_Solution/Activities/Assignment_Solution_TMB_plot.pdf"
        File SignatureCount = "~{OutputFolder}/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
    }

    runtime {
        docker: "us.gcr.io/tag-public/sigprofiler:v1"
        memory: "8 GB"
        disks: "local-disk 20 HDD"
    }
}


task PlotSignatures {
  input {
    File SignatureCount
  }

  command {
    python3 <<EOF
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    SigCounts = pd.read_csv("${SignatureCount}", sep='\t', header=0)
    SigCounts = pd.DataFrame(SigCounts)

    # Calculate proportions
    signature_cols = SigCounts.columns[1:]  # Exclude the 'Samples' column
    SigCounts[signature_cols] = SigCounts[signature_cols].div(SigCounts[signature_cols].sum(axis=1), axis=0)

    # Reshape the data
    SigCounts_long = SigCounts.melt(id_vars=["Samples"], var_name="Signature", value_name="Proportion")
    SigCounts_long = SigCounts_long[SigCounts_long["Proportion"] > 0]

    # Plot the data
    plt.figure(figsize=(16, 9))
    sns.scatterplot(data=SigCounts_long, x="Samples", y="Signature", size="Proportion", sizes=(20, 200), legend=False)
    plt.xticks(rotation=90)
    plt.xlabel("Sample Name", fontsize=16)
    plt.ylabel("Signature", fontsize=16)
    plt.title("Signature Proportions by Sample", fontsize=20, pad = 20)
    plt.grid(axis='y')
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    plt.tight_layout()
    plt.savefig("signature_proportions.pdf", format="pdf")
    EOF
  }

  output {
    File signature_proportions_pdf = "signature_proportions.pdf"
  }

  runtime {
        docker: "us.gcr.io/tag-public/sigprofiler:v1"
        memory: "8 GB"
        disks: "local-disk 20 HDD"
  }
}
