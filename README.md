# wgbs-genetic-variants
Instructions for genetic variant detection from whole genome bisulfite sequencing (WGBS) data

## Introduction

## Results of Our Accuracy Analysis

## Instructions for Running Revelio
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
## WGS variant calling with Sarek using Freebayes
1)	Sequencing quality control and trimming of raw fastq files using FASTQC and FASTP
2)	Map reads to reference using BWAMEM1
3)	Process BAM file (GATKMarkduplicates)
4)	Call variants using the Freebayes variant caller (SNPs and Indels)
```
nextflow run /home/.nextflow/assets/nf-core/sarek/main.nf  --input /home/sarek/in4276_freebayes_wgs//home/hpatel/sarek/in4319/in4319_input.csv -profile awsbatch --awsregion us-east-1 --awsqueue 'arn:aws:batch:us-east-1:002226384833:job-queue/methylseq'  --step variant_calling --outdir 's3://zymo-filesystem/home/hpatel/sarek/in4319_wgs_freebayes/' -w 's3://zymo-filesystem/tmp/hpatel/in4319_work/' --genome null --igenomes_ignore --fasta 's3://zymo-filesystem/home/hpatel/reference_genomes/GIAB/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta' --skip_tools baserecalibrator --tools freebayes
```
5)	Filtering VCF file to retain high-quality variants. QUAL > 30 and DP > 10 seem to be optimal for filtering.

`bcftools filter -i 'QUAL > 30 && FORMAT/DP > 10' -O z -o "$output_dir/filtered1.vcf.gz" "$input_file1"`



## Instructions for Accuracy (Precision/Recall) Analysis of GIAB Sample

**GIAB HG38 Fasta:**

Make the fai file with `samtools faidx`

[https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38\_GIABv3\_no\_alt\_analysis\_set\_maskedGRC\_decoys\_MAP2K3\_KMT2C\_KCNJ18.fasta.gz](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz "‌")

**precisionFDA truth dataset for HG001 (VCF + stratification BED):**

[https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv3.3.2/GRCh38/](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv3.3.2/GRCh38/ "smartCard-inline")

[**rep.py**](http://rep.py "‌")**:**

Must run with python2 + jinja2.

[https://github.com/ga4gh/benchmarking-tools/blob/master/reporting/basic/bin/rep.py](https://github.com/ga4gh/benchmarking-tools/blob/master/reporting/basic/bin/rep.py "smartCard-inline")

sample command:

```
python /home/mjin/scratch/projects/GIAB_Fixture/benchmarking-tools/reporting/basic/bin/rep.py \
  in4307-1_hap.py:../happy/in4307_1.roc.all.csv.gz \
  in4307-2_hap.py:../happy/in4307_2.roc.all.csv.gz \
  in4307-3_hap.py:../happy/in4307_3.roc.all.csv.gz \
  in4307-4_hap.py:../happy/in4307_4.roc.all.csv.gz \
  in4307-5_hap.py:../happy/in4307_5.roc.all.csv.gz \
  in4307-6_hap.py:../happy/in4307_6.roc.all.csv.gz \
  in4307-7_hap.py:../happy/in4307_7.roc.all.csv.gz \
  in4307-8_hap.py:../happy/in4307_8.roc.all.csv.gz \
  -o comparisons.html
```

[**hap.py**](http://hap.py "‌")**:**

Don’t install the program from bioconda. It won’t work because of incompatible C compiler. Use the docker container:

`docker run -it --rm -v <host>:<container> quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0 bash`

sample command (inside container):

```
hap.py \
  /home/mjin/scratch/projects/GIAB_Fixture/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
  ../variant/in4307_1.deepvariant.vcf.gz \
  -f /home/mjin/scratch/projects/GIAB_Fixture/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed \
  -r /home/mjin/scratch/projects/GIAB_Fixture/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta \
  -o in4307_1
```
