# Variant Calling from Whole Genome Bisulfite Sequencing Data (WGBS): From FASTQ to VCF

This tutorial provides a step-by-step guide on how to perform variant calling on Whole Genome Bisulfite Sequencing (WGBS) data, from BAM files to obtaining a VCF file. We will use the following tools:

1. **Revelio** to convert WGBS BAM files to genomics-compatible BAM files.
2. **FreeBayes** for variant calling.

For the purposes of this tutorial, we already assume your WGBS data has been appropriately trimmed, aligned (for example, using a tool such as Bismark), and QC-checked. We will begin with an aligned BAM file and demonstrate how to obtain called variants in a VCF (Variant Call Format) file.

## Conceptual Introduction

We have elected to use a tool called Revelio, which transforms a BAM file from WGBS data into one that is compatible with the same tools as any genomics analysis. The advantage is that standard variant calling tools can now be used on WGBS data. Also, if a researcher is performing both WGBS and whole genome sequencing (WGS, or DNA-seq), the same variant calling pipeline could be used to compare between the datasets. Revelio achieves this by transforming quality scores based on the likelihood that a C->T conversion is due to a genetic change or a bisulfite conversion. It is able to do this, for stranded (directional), WGBS libraries by leveraging unconverted information present on the second strand. The quality score conversion will effectively mask any bisulfite conversion events, but will allow the variant caller to analyze any other mutation events. See the Revelio publication for more information on how the algorithm works.

## Step 1: Running Revelio
**refer to [Revelio's repository](https://github.com/bio15anu/revelio) for installation instructions**

**generate MD tags based on sample bam and reference genome fasta**
```
samtools calmd -b sample_1.bam genome.fa 1> sample_1_calmd.bam 2> /dev/null
samtools index calmd.bam
```
**generate masked bam**
```
./revelio.py sample_1_calmd.bam sample_1_masked.bam
samtools index sample_1_masked.bam
```
## Step 2: WGS variant calling with Sarek using Freebayes
1)	Sequencing quality control and trimming of raw fastq files using FASTQC and FASTP
2)	Map reads to reference using BWAMEM1
3)	Process BAM file (GATKMarkduplicates)
4)	Call variants using the Freebayes variant caller (SNPs and Indels)
```
nextflow run /home/.nextflow/assets/nf-core/sarek/main.nf  \
--input ./input.csv \
-profile awsbatch \
--awsregion us-east-1 \
--awsqueue 'methylseq'  \
--step variant_calling \
--outdir './sarek/wgs_freebayes/' \
-w './work/' \
--genome null \
--igenomes_ignore \
--fasta './reference_genomes/GIAB/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta' \
--skip_tools baserecalibrator \
--tools freebayes
```
5)	Filtering VCF file to retain high-quality variants. QUAL > 30 and DP > 10 seem to be optimal for filtering.

`bcftools filter -i 'QUAL > 30 && FORMAT/DP > 10' -O z -o "$output_dir/filtered1.vcf.gz" "$input_file1"`

## Notes

- Ensure all tools are properly installed and accessible in your PATH.
- Adjust file paths and parameters as necessary for your specific data and computational environment.
- The quality of the input data and the parameters used for alignment and variant calling can significantly affect the results.