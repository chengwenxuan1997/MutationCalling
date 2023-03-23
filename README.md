# Call somatic mutation from matched tumor and normal sample (long reads WES sample)
Author: Wenxuan Cheng (chengwenxuan1997@outlook.com)  
Date: 2023-03-21

# 1: Base usage for the program (for just using the pipeline)
## for a new project
The core script is only "script.sh"
module: bedtools, picard, gatk4, perl
Only the InputData/sample_sheet.csv need to be modifed (And bam file should be put under InputData/)
	InputData/sample_sheet.csv looks like:

	#TumorBam,NormalBam
	patient2.tumor.bam,patient2.normal.bam
if the line starts with #, this pair of samples will be discarded.

## Required reference file list
Genome: Reference/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa  
dbsnp: Reference/GRCh37.p13.dbSNP151.All.vcf.gz  
	download from https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz  
VEP Cache: Reference/homo_sapiens_vep_109_GRCh37  
	download from https://ftp.ensembl.org/pub/grch37/release-109/variation/indexed_vep_cache/homo_sapiens_vep_109_GRCh37.tar.gz  
blacklist: Reference/hg19-blacklist.v2.bed  
	download from https://github.com/Boyle-Lab/Blacklist/tree/master/lists  
	hg19, hg38, mm10 are all provided under Reference/  
commonVCF: Reference/GetPileupSummaries/small_exac_common_3_b37.vcf.gz  
	download from ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2  
	hg19(b37), hg38 are all provided under Reference/  
	**there is no available file for the mouse species, so it should be "null" if it is a mouse sample**  
gnomad: Reference/af-only-gnomad.raw.sites.b37.vcf.gz  
	download from ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2  
	hg19(b37), hg38 are all provided under Reference/GetPileupSummaries  
	**there is no available file for the mouse species, so it should be "null" if it is a mouse sample**  

## Output list (Under the Results/)  
Merged.vcf: the merged vcf for all sample  
Merged.maf: the merged maf file required for downstream analysis  
*.final.vcf: the mutation calling output for each sample  
*.final.vcf: the annotated vcf file for each sample  
*.anno.maf: the maf file for each sample  


# 2: More about the program (for understanding the whole pipeline)  
The main body program consists of three parts: ProcessBam.sh, SomaticVariantCalling.sh and VariantAnnotation.sh  
After preprocessing all bam files, the script will call mutations with MuSE, MUTECT2, VARSCAN, and then intersect the three mutation sets. Finally, the script will annotate the all-sample merged vcf with VEP and output a maf file for downstream analysis.  

## 2.1: ProcessBam.sh (Preprocess the bam file)  

### Description  
The script performs Base Quality Score Recalibration (BQSR) for each bam file and somatic variant calling for each normal bam file (for building a normal panel).   
Because there is no "chr" prefix in the common variant database (including dbsnp, gnomad, gatk), the ProcessBam.sh will provide two ways to clean the bam (if there is a "chr" prefix in the bam). One is to remove them with "sed 's/chr//g'" simply if the realignment is set as "F". The other is to realign these reads to provided GRCh genome if the realignment is set as "T".   
Besides the script does not help in handling discordance of chromosome size between alignment genome and variant database if they were built with different versions (even minor version differences, GRCH37.p13, and GRCH37.p10, rather than GRCH37 and GRCh38). Thus it is highly possible to report a error, and the best solution is to realign these bam to newest reference genome or to switch the matched variant database (it is difficult to find suitable version for the database from gnomad and gatk).  

### Argument list
-b, bam: the file path of the bam to be process  
-r, Reference: the file path of the reference genome, end up with "fa"  
-t, thread: the number of thread for alignment, default is 10  
-i, id: the sample identifier for this bam  
-d, dbsnp: the file path of dbsnp vcf  
-a, realignment: whether realign to genome without "chr" prefix. If "T" perform the realignment, and if "F" remove the "chr" prefix with sed command  
-c, class: whether the sample belongs to tumor or normal. If "Normal" variant calling will be performed for a panel of normals (PoN), otherwise nothing will be done.  

## 2.2: SomaticVariantCalling.sh (call the mutation from paired bam files)

### Description
The script will call mutations with MuSE, MUTECT2 and VARSCAN for each paired sample (strelka written with python2 is excluded for its conflicts with picard although it is a powerful tool).  

### Argument list
-t, TumorBam, the file path of tumor bam  
-n, NormalBam, the file path of normal bam  
-r, Reference, the file path of Reference genome  
-d, dbsnp, the file path of dbsnp vcf  
-i, ID, the identifier for this paired of samples  
-g, gnomad: list for comman variants  
-c, common VCF: list for comman variants  

## 2.3: VariantAnnotation.sh (annotate the merged vcf file)

### The description 
The script will merge vcfs from all samples, annotate the merged vcf with VEP (Variant Effect Predictor) and turn the annotated vcf into maf-format.  

### Argument list
-c, CachePath, the file path of VEP cache  
-r, Reference, the file path of genome  

