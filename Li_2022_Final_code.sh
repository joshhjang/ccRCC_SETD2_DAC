#Software versions
hisat2-2.1.0
trim_galore-0.6.0
subread-2.0.0
picard-2.23.3
python-3.8.1
Stringtie-2.1.3b
ggsashimi-1.0.0
gffcompare-0.12.2
R-4.0.2
bwa-0.7.17
samtools-1.9
macs2-2.2.7.1
HOMER-4.11.1
idr-2.0.3
deeptools-3.4.3

##############     Processing totalRNA-seq reads ###################
#Trim SE reads
trim_galore --fastqc -o ../trimmedreads $Sample_fastq_R1 > $Sample_trimgalore_summary.txt 
#Align totalRNA-seq SE reads with HISAT2
hisat2 -x $HISAT2_reference --rna-strandness R --dta -U $Sample_fastq_R1_trimmed -p 16 --summary-file $Sample_HISAT_summary.txt | samtools view -bS - | samtools sort -T $Sample.tmp -O bam -o $Sample.bam 

###Trim PE reads
trim_galore --fastqc -o ../trimmedreads --paired $Sample_fastq_R1 $Sample_fastq_R2 > $Sample_trimgalore_summary.txt 
#Align totalRNA-seq PE reads with HISAT2
hisat2 -x $HISAT2_reference --rna-strandness RF --dta -1 $Sample_fastq_R1_trimmed -2 $Sample_fastq_R2_trimmed -p 16 --summary-file $Sample_HISAT_summary.txt | samtools view -bS - | samtools sort -T $Sample.tmp -O bam -o $Sample.bam

#Index bam files 
samtools index $Sample.bam

#Quantify TE expression with FeatureCounts (SE reads)
featureCounts -s 2 -T 14 -t exon -g transcript_id -f -a GRCh38_GENCODE_rmsk_TE.gtf -o $Sample_rpms.tsv $Sample.bam --primary -Q 20 --ignoreDup
#Quantify TE expression with FeatureCounts (PE reads)
featureCounts -s 2 -T 16 -t exon -g transcript_id -p -f -a GRCh38_GENCODE_rmsk_TE.gtf -o $Sample_rpms.tsv $Sample.bam --primary -Q 20 --ignoreDup

#Quantify Genecode V37 gene expression with stringtie
stringtie $Sample.bam -p 14 -o $Sample.gtf -G /gencode.v37.annotation.gtf -A $Sample_gene_abundance.txt -e -b $Sample_ballgown

#Make gene counts file using stringtie prepDE3.py
python3 prepDE3.py -i sample_lst.txt

#Perform De Novo transcript assembly and merge with Stringtie and annotate transcript isoforms with gffcompare
stringtie $Sample.bam -p 16 -l $Sample -o $Sample_denovoAssembly.gtf --fr
stringtie --merge -c 2 -f 0.05 -i -p 16 --rf -o $Sample_denovo_merged.gtf -G gencode.v37.annotation.gtf $Sample_assembly_GTF_list.txt
gffcompare -r gencode.v37.annotation.gtf -o $Sample_gffcompare_V37_denovo $Sample_denovo_merged.gtf



############## Processing ChIP-seq reads ################
###Trim PE reads
trim_galore --fastqc -o ../trimmedreads --paired $Sample_fastq_R1 $Sample_fastq_R2 > $Sample_trimgalore_summary.txt 
#Align ChIP-seq PE reads with HISAT2
bwa mem -t 8 hg38.fa $Sample_fastq_R1_trimmed -2 $Sample_fastq_R2_trimmed | samtools sort -@8 -o $Sample.sorted.bam -
#Remove duplicate reads (picard-2.23.3)
java -Xms8g -Xmx16g -Djava.io.tmpdir=./tmp -jar $PICARD MarkDuplicates I=$Sample.sorted.bam O=$Sample.sorted.rmdup.bam M=$Sample_dup_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true

## Call histone peaks using MACS2
macs2 callpeak -t $Sample.sorted.rmdup.bam -c Control_Input_combined.sorted.rmdup.bam -n $Sample -g 2.7e9 --outdir ./ -f BAMPE -s 50 --broad --broad-cutoff 0.1 #Broad peaks
macs2 callpeak -t $Sample.sorted.rmdup.bam -c Control_Input_combined.sorted.rmdup.bam -n $Sample -g 2.7e9 --outdir ./ -f BAMPE -s 50 #Narrow peaks

#Call IDR peaks
idr --samples $Sample_Rep1_peaks.Peak $Sample_Rep2_peaks.Peak --input-file-type "PeakType" -o $Sample.IDR_narrowPeaks.txt --log-output-file $Sample_narrowPeaks.IDR.log --idr-threshold 0.05 --plot --verbose

#Annotate genomic region of peaks with Homer
annotatePeaks.pl $Sample_peak.bed hg38 -genomeOntology $Sample_genomeOntology -gtf gencode.v37.annotation.gtf > $Sample_annotatePeaks.bed

#Generate Average histone signals across Gene Body
computeMatrix scale-regions --missingDataAsZero -bs 100 -p 12 -R gencode.v37.AutosomalTranscriptsOnly.forDeepTools.gtf -S $Sample_EV_Histone_Day0.bw $Sample_EV_Histone_Day5.bw $Sample_EV_Histone_Day23.bw $Sample_KO_Histone_Day0.bw $Sample_KO_Histone_Day5.bw $Sample_KO_Histone_Day23.bw -m 5000 -b 3000 -a 3000 -out $Histone_GeneBody_O786_combined.tab.gz 
plotHeatmap  --heatmapWidth 6 -m $Histone_GeneBody_O786_combined.tab.gz --perGroup -out $Histone_GeneBody_O786_combined.pdf --missingDataColor "grey" --colorMap "Color" --heatmapHeight 10

















