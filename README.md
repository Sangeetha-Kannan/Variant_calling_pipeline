# Genomic Data Processing Pipeline
Welcome to the comprehensive pipeline for processing and analyzing high-throughput sequencing data! This workflow takes you from raw sequencing reads all the way to quality-controlled, ready-for-variant-calling BAM files. Leveraging powerful tools like BWA, SAMtools, and GATK, this pipeline ensures accurate and efficient data processing, making it a robust solution for any genomic analysis project.

Here's how you can perform each step from creating directories, downloading data, performing quality checks and proceeding with analysis.

These steps provide a comprehensive approach to germline variant calling, genotyping, and annotation, essential for genetic research and understanding heritable genetic variations.

####1. Setting up Directories
##### Step 1: Create a Directory for the Reference Genome 
Create a directory such as ref_genome where you will store the reference genome files. Navigate into this directory.
```shell
$ mkdir ref_genome
$ cd ref_genome
```
####2. Download and Decompress Reference Genome
##### Step 2: Download Reference Genome from UCSC genome Browser
Download the reference genome file (here, Chromosome 6) from UCSC Genome browser and decompress it using **gunzip**. 
```shell
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr6.fa.gz
gunzip chr6.fa.gz
```
####3. Configure SRA Toolkit
#####Step 3: Download and Configure SRA Toolkit
Navigate to the home directory or any preferred location outside ref_genome and download the SRA Toolkit and configure it by extracting the files.  
Navigate into the bin directory of  SRA Toolkit and execute ./vdb-config to open the configuration tool. Here, set up the environment and permissions for using SRA Toolkit (go to TOOLS option and set to current directory)
```shell
   cd ~
  wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
  tar -xzvf sratoolkit.current-ubuntu64.tar.gz
  sratoolkit.3.1.0-ubuntu64 sra_toolkit
  cd sra_toolkit/bin
  ./vdb-config -i
```
####4. Set up the Raw Data Directory
#####Step 4: Create a Directory for Raw Data such as raw_data where you will store the FASTQ files. Navigate into this directory
```shell
mkdir raw_data
cd raw_data
```
####5. Download FASTQ files
#### Step5: Download the paired-end FASTQ files from the SRA database
Runs the fastq-dump command from the SRA Toolkit to download and split the paired-end reads ERR11468775 into separate files and also compress the output FASTQ files using gzip

```shell
./sra_toolkit/bin/fastq-dump --split-files --gzip ERR11468775 -X 50000
```
####6. Install and Run FASTQC
##### Step 6: Installs FastQC, a quality control tool for high throughput sequence data
Note that you may need to additionally add --fix-missing option during installation
```shell
sudo apt-get update
sudo apt-get -y install fastqc
fastqc ERR11468775_1.fastq.gz ERR11468775_2.fastq.gz
```
####7. Install and Run FASTP
##### Step 7: Install FASTP and use it to trim and filter the raw FASTQ files
```shell
wget http://opengene.org/fastp/fastp
chmod a+x fastp
./fastp -i ERR11468775_1.fastq.gz -o ERR11468775_trimmed_1.fastq.gz -I ERR11468775_2.fastq.gz -O ERR11468775_trimmed_2.fastq.gz --detect_adapter_for_pe -f 10 -g -l 50 -c -h ERR11468775_fastp.html -w 10
```
####8. Install BWA 
##### Step 8: Installation of BWA for Read Mapping
Navigate to the home repository. Install BWA, a software package for mapping low-divergent sequences against a large reference genome. It is process of aligning reads to the reference genome.
```shell
git clone https://github.com/lh3/bwa.git
cd bwa; make
```
##### Step 9: Indexing the Reference genome
Indexing Reference genome before mapping reads help create auxillary data structures that allow efficient access and rapid alignment of reads.
```shell
./bwa index -a bwtsw -p ../ref_genome/chr6_ref ../ref_genome/chr6.fa
```
##### Step 10: Read Mapping
The sequencing reads are compared to the indexed reference genome, and the best matching positions are identified for each read.
```shell
bwa mem ../ref_genome/chr6_ref ../raw_data/ERR11468775_trimmed_1.fastq.gz ../raw_data/ERR11468775_trimmed_2.fastq.gz -t 10 -o sample1.sam
```
####9.  Conversion of SAM to BAM file
Convert the SAM file (text format) to BAM format (binary format), which is more efficient for storage and processing.
##### Step 11: Ensure files are in the right directory
Your reference genome chr6.fa is in ref_genome directory and sample1.sam is in bwa directory. First, we need to copy these files into your current working directory. Create a new directory here or navigate to the home directory.
```shell
cp ../ref_genome/chr6.fa .
cp ../bwa/sample1.sam .
```
##### Step 12: Install Samtools
Samtools is a suite of programs for interacting with high-throughput sequencing data that enables conversion and manipulation of SAM/BAM files
```shell
sudo apt-get update
sudo apt-get -y install samtools
```
##### Step 13: Convert SAM to BAM
```shell
samtools view -bo sample1.bam sample1.sam
```
####10. Install docker to run GATK container
(If you are using Gitpod, it has docker installed by default)
##### Step 14: Download GATK Docker Image
Pull the GATK Docker image from the Broad Institute's repository.
```shell
sudo docker pull broadinstitute/gatk:latest
```

##### Step 15: Download Known Variant files
Download the VCF files containing known genetic variants and their associated index files.
```shell
Wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/genomics-public-data/resources//hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
```
##### Step 16: Run Docker in Interactive Mode
Run the GATK Docker container in interactive mode, mounting the current directory to /data inside the container for easy access to your files
```shell
sudo docker run -it -v $PWD:/data broadinstitute/gatk:latest
```
##### Step 17: Add Read Group Information
Add read group information to the BAM file. Read groups help differentiate between different sequencing runs or samples
Add metadata to the BAM file that identifies the read group, library, platform, platform unit, and sample. This is crucial for downstream analyses.
```shell
gatk AddOrReplaceReadGroups -I sample1.bam -O sample1_withRG.bam -ID 1 -LB lib1 -PL ILLUMINA -PU unit1 -SM sample1
```
##### Step 18: Sort BAM file
Sort the BAM file by coordinates.
Ensure the BAM file is sorted by genomic coordinates, which is required for many downstream processes like variant calling and visualization.
```shell
gatk SortSam -I sample1_withRG.bam -O sorted_sample1.bam -SO coordinate
```
##### Step 19: Check Mapping Results
Generate statistics about the alignments in the BAM file. Get a summary of alignment statistics, such as the total number of reads, mapped reads, and duplicate reads, to assess the quality of the alignment.
```shell
samtools flagstat sorted_sample1.bam
```
##### Step 20: Mark Duplicates
Identify and mark duplicate reads in the BAM file.
Mark duplicate reads, which are often artifacts from the PCR amplification process. This helps to prevent these duplicates from biasing downstream analyses.
```shell
gatk MarkDuplicates -I sorted_sample1.bam -O markedDups.bam -M metrics_duplicates
```
##### Step 21: Base Quality Score Recalibration (BQSR)
Perform BQSR to correct systematic errors made by the sequencer when estimating the quality score of each base call.
###### a) Create a sequence dictionary for reference genome
This command generates a sequence dictionary for the reference genome chr6.fa. The dictionary file is necessary for various GATK tools to quickly access the sequences in the reference genome
```shell
gatk CreateSequenceDictionary -R chr6.fa
```
###### b)	Index the Reference Genome with Samtools
This command indexes the reference genome chr6.fa using samtools. The index file allows for efficient access to the genomic positions in the reference genome
```shell
samtools faidx chr6.fa
```
###### c)	Generate a recalibration table based on known variants
This command generates a recalibration table (sample1_recal_data.table) for the input BAM file (markedDups.bam) using the reference genome (chr6.fasta) and known variant sites (Homo_sapiens_assembly38.dbsnp138.vcf). The recalibration table contains information about systematic errors in the base quality scores of the input BAM file
```shell
gatk BaseRecalibrator -I markedDups.bam -R chr6.fasta --known-sites Homo_sapiens_assembly38.dbsnp138.vcf -O sample1_recal_data.table
```
#####d) Apply BQSR to bam file
This command applies the base quality score recalibration to the input BAM file (sample1_markedDups.bam) using the recalibration table (sample1_recal_data.table). The output is a recalibrated BAM file (sample1_recal.bam), which has adjusted base quality scores that correct for systematic errors
```shell
../../gatk/gatk ApplyBQSR -R ../ref_genome/chr6.fa -I sample1_markedDups.bam --bqsr-recal-file ../sample1_recal_data.table -O sample1_recal.bam
```
##### Step 22: Indexing Recalibrated BAM File
This command creates an index file for the recalibrated BAM file (sample1_recal.bam). Indexing is essential as it allows for faster access and retrieval of data from the BAM file during downstream analyses, such as variant calling
```shell
samtools index /data/mapping/sample1_recal.bam
```
####11. Calling Somatic Variants Using GATK
#####a) Running Tumor-Only Pipeline
This command uses GATK's Mutect2 to perform somatic variant calling. It analyzes the recalibrated BAM file (sample1_recal.bam) against the reference genome (chr6.fa) to identify somatic mutations. The result is a VCF file (sample1.vcf.gz) that contains the list of potential somatic variants.
```shell
gatk Mutect2 -I sample1_recal.bam -R ../ref_genome/chr6.fa -O sample1.vcf.gz
```
#####b) Filtering Variants
This command filters the called somatic variants using GATK's FilterMutectCalls. It evaluates the variants in the VCF file (sample1.vcf.gz) against the same reference genome to apply various filters that help in removing false positives and refining the list of somatic variants. The filtered VCF file (filtered_sample1.vcf.gz) contains the high-confidence somatic variants after applying the necessary filters.
```shell
gatk FilterMutectCalls -R ../ref_genome/chr6.fa -V sample1.vcf.gz -O filtered_sample1.vcf.gz
```
##### c) Extracting Indels
This command selects only insertion and deletion mutations (indels) from the VCF file (sample1.vcf.gz). The result is appended to sample1_indels.vcf, a file specifically containing indel variants from the sample.
```shell
bcftools view --types indels sample1.vcf.gz >> sample1_indels.vcf
```
##### d) Extracting SNPs
Similar to the previous command but tailored for single nucleotide polymorphisms (SNPs). It extracts only SNP-type variants from the VCF file.
```shell
bcftools view --types snps sample1.vcf.gz >> sample1_snps.vcf
```
##### e) Filtering VCF File for Quality
This command filters the VCF file (sample1.vcf.gz) to include only variants with a quality score greater than 50. The quality score (%QUAL) reflects the confidence in the correctness of a variant call.
```shell
bcftools filter -i '%QUAL>50' sample1.vcf.gz
```
### 12. Alternatively, Germline Variant Calling
#####Step 1: Download Reference Files
These commands download the HapMap reference files in VCF format along with their index files. HapMap data are used as a reference in variant calling processes to help distinguish between common population variants and potential novel or deleterious variants
```shell
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi
```
##### Step 2: Variant Calling with HaplotypeCaller
These commands run the GATK HaplotypeCaller to perform germline variant calling on each sample (sample1: father, sample2: mother, sample3: child). Each BAM file is processed to detect variants that are then stored in a gVCF file for each sample.
```shell
gatk HaplotypeCaller -R chr6.fa -I sample1_recal.bam -O sample1_germline.g.vcf.gz -ERC GVCF
gatk HaplotypeCaller -R chr6.fa -I sample2_recal.bam -O sample2_germline.g.vcf.gz -ERC GVCF
gatk HaplotypeCaller -R chr6.fa -I sample3_recal.bam -O sample3_germline.g.vcf.gz -ERC GVCF
```
##### Step3: Importing VCF to GenomicsDB
This command integrates the gVCF files into a GenomicsDB datastore. This integration is crucial for managing and querying large genomic datasets efficiently, particularly when performing joint genotyping on multiple samples.
```shell
gatk GenomicsDBImport -V sample1_germline.g.vcf.gz -V sample2_germline.g.vcf.gz -V sample3_germline.g.vcf.gz --genomicsdb-workspace-path test_db --intervals chr6
```
##### Step 4: Genotyping
This command performs joint genotyping on the variant data stored in the GenomicsDB workspace. It aggregates the gVCF data from all samples and calculates the most likely genotypes.
```shell
gatk GenotypeGVCFs -R chr6.fa -V gendb://test_db -O try_new.vcf.gz
```
### Annotation of Variants with VEP (Variant Effect Predictor)
 For annotating the variants, you can use the Variant Effect Predictor (VEP) tool provided by ENSEMBL. VEP annotates variants from VCF files to explore their potential biological significance.
