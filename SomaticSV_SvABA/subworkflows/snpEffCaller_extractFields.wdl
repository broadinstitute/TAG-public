#snpEffCaller
#snapshot 13
#snpEff for annotation of vcf files using hg19. 

workflow SnpEff {
		File snpEff_jar
        File snpSift_jar
        File inputVCF
        File configFile
        File databaseFile
        String outputName
        
        call snpEffCaller {
         input:
           snpEff_jar = snpEff_jar,
           snpSift_jar = snpSift_jar,
           inputVCF = inputVCF,
           configFile = configFile,
           databaseFile = databaseFile,
           outputName = outputName,
        }
}

task snpEffCaller {
        File snpEff_jar
        File snpSift_jar
        File inputVCF
        File databaseFile
        String outputName
        File configFile

        command {
        		subdir=$(dirname ${snpEff_jar})
        		unzip ${databaseFile} -d $subdir
                java -Xmx3g -jar ${snpEff_jar} -v -config ${configFile} -canon -nodownload -hgvs1LetterAa hg19 ${inputVCF} > ${outputName}_snpEff.vcf
                
                java -Xmx3g -jar ${snpSift_jar} extractFields -e "." ${outputName}_snpEff.vcf CHROM POS ID REF ALT QUAL FILTER SCTG "ANN[*].GENE" "ANN[*].GENEID" > ${outputName}_snpEff_snpSift.vcf
        }

        output {
                File vcfOut = "${outputName}_snpEff.vcf"
                File vcfOutExtractFields = "${outputName}_snpEff_snpSift.vcf"
        }
        
		runtime {
        	docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
            memory: "3.75GB"
            cpu: 1
            disks: "local-disk 200 HDD"
    }
}