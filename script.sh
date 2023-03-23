# 
# included module: bedtools, picard, gatk4, perl #(removed: MuSE, VarScan, Strelka2)
# gatk version: 4.3.0 (some arguments differ between different version of gatk)
# about the
# use bam without prefix "chr"


Reference="Reference/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa"    # reference genome
dbsnp="Reference/GRCh37.p13.dbSNP151.All.vcf.gz"                        # dbsnp vcf
CachePath="Reference/homo_sapiens_vep_109_GRCh37"                       # VEP cache
blacklist="Reference/hg19-blacklist.v2.bed"                             # black list
gnomad="Reference/af-only-gnomad.raw.sites.b37.vcf.gz"                  # gnomad
commonVCF="Reference/GetPileupSummaries/small_exac_common_3_b37.vcf.gz" # gatk common vcf

# create index for genome
## Create index for reference genome
cmd="bwa index $Reference 2>/dev/null"
if [ ! -f $Reference.ann ]; then echo $cmd; eval $cmd; fi
cmd="samtools faidx $Reference 2>/dev/null"
if [ ! -f $Reference.fai ]; then echo $cmd; eval $cmd; fi
cmd="picard CreateSequenceDictionary R=$Reference O=$(echo $Reference|sed 's/fa/dict/g') 2>/dev/null"
if [ ! -f $(echo $Reference|sed 's/fa/dict/g') ]; then echo $cmd; eval $cmd; fi

# create index for dbsnp vcf
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz
if [ ! -f "$dbsnp.tbi" ];then
  echo "$dbsnp.tbi not exists"
  cmd="gatk IndexFeatureFile --input $dbsnp --QUIET true 2>/dev/null"
  echo $cmd; eval $cmd
fi

## vep cache
# https://ftp.ensembl.org/pub/release-109/variation/indexed_vep_cache/homo_sapiens_vep_109_GRCh38.tar.gz
# https://ftp.ensembl.org/pub/grch37/release-109/variation/indexed_vep_cache/homo_sapiens_vep_109_GRCh37.tar.gz

# create interval_list (whole genome/exome subtract the blacklist)
awk -v FS="\t" -v OFS="\t" '{print $1 FS "0" FS ($2-1)}' $Reference.fai > $Reference.bed
bedtools subtract -a $Reference.bed -b $blacklist > $Reference.interval.bed
cmd="gatk BedToIntervalList -I $Reference.interval.bed -O $(echo $Reference|sed 's/fa/interval_list/g')"
cmd="$cmd -SD $(echo $Reference|sed 's/fa/dict/g') 2>/dev/null"
echo $cmd; eval $cmd


echo "=============================================="
echo "=========Step 1: process the bam file========="
echo "=============================================="

for i in $(cat InputData/sample_sheet.csv|grep -v \#|wc -l|xargs seq)
do
{
  tmp=$(cat InputData/sample_sheet.csv|grep -v \#|sed -n ${i}p)
  TumorBam="InputData/$(echo $tmp|cut -d, -f1)"       # Tumor file
  NormalBam="InputData/$(echo $tmp|cut -d, -f2)"      # Normal file
  TumorID=$(samtools view -H $TumorBam|grep ^@RG|sed "s/\t/\n/g"|grep SM|sed 's/SM://g')
  NormalID=$(samtools view -H $NormalBam|grep ^@RG|sed "s/\t/\n/g"|grep SM|sed 's/SM://g')

  echo
  echo "==============$TumorID==$NormalID==========="

  cmd1="Codes/1.ProcessBam.sh -b $TumorBam -r $Reference -i $TumorID -d $dbsnp -t 10 -a F"
  cmd2="Codes/1.ProcessBam.sh -b $NormalBam -r $Reference -i $NormalID -d $dbsnp -t 10 -a F -c Normal"

  echo $cmd1; eval $cmd1 &
  echo $cmd2; eval $cmd2 &
  wait

  echo "============================================"
  echo
  echo
}
done

echo "==============================================="
echo "=======Step 2: Call the somatic Mutation======="
echo "==============================================="

cmd="gatk GenomicsDBImport -R $Reference -L $(echo $Reference|sed 's/fa/interval_list/g')"
cmd="$cmd  --genomicsdb-workspace-path pon_db -V $(ls|grep normal.vcf.gz$|sed ':a;N;$!ba; s/\n/ -V /g') 2>/dev/null"
echo $cmd; eval $cmd

cmd="gatk CreateSomaticPanelOfNormals -R $Reference --germline-resource $gnomad -V gendb://pon_db -O pon.vcf.gz 2>/dev/null"
echo $cmd; eval $cmd

for i in $(cat InputData/sample_sheet.csv|grep -v \#|wc -l|xargs seq)
do
{
  tmp=$(cat InputData/sample_sheet.csv|grep -v \#|sed -n ${i}p)
  TumorBam="InputData/$(echo $tmp|cut -d, -f1)"       # Tumor file
  NormalBam="InputData/$(echo $tmp|cut -d, -f2)"      # Normal file
  TumorID=$(samtools view -H $TumorBam|grep ^@RG|sed "s/\t/\n/g"|grep SM|sed 's/SM://g')
  NormalID=$(samtools view -H $NormalBam|grep ^@RG|sed "s/\t/\n/g"|grep SM|sed 's/SM://g')
  
  echo
  echo "===========$TumorID==$NormalID========="

  cmd="Codes/2.SomaticVariantCalling.sh -t $TumorID.clean.calib.bam -n $NormalID.clean.calib.bam"
  cmd="$cmd -r $Reference -d $dbsnp -i $TumorID -g $gnomad -c $commonVCF"
  echo $cmd; eval $cmd
  
  ls|grep $TumorID|grep -v .final.vcf|xargs rm
  ls|grep $NormalID|grep -v .final.vcf|xargs rm

  echo "======================================="
  echo
  echo
}
done

wait

echo "==============================================="
echo "=========Step 3: annotate the mutation========="
echo "==============================================="

cmd="Codes/3.VariantAnnotation.sh -c $CachePath -r $Reference"
echo $cmd; eval $cmd

rm -r pon_db
ls|grep pon|xargs rm