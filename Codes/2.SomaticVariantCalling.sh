# Somatic Variant Calling

TumorBam="patient2_T.clean.calib.bam"
NormalBam="patient2_N.clean.calib.bam"
Reference="Reference/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa"
dbsnp="Reference/GRCh37.p13.dbSNP151.All.vcf.gz"
ID="patient2"
gnomad="null"
commonVCF="null"

while getopts t:n:r:d:i:g:c: flag
do
    case "${flag}" in
        t) TumorBam=${OPTARG};;
        n) NormalBam=${OPTARG};;
        r) Reference=${OPTARG};;
        d) dbsnp=${OPTARG};;
        i) ID=${OPTARG};;
        g) gnomad=${OPTARG};;
        c) commonVCF=${OPTARG}
    esac
done
# echo "TumorBam: $TumorBam";
# echo "NormalBam: $NormalBam";
# echo "Reference: $Reference";
# echo "dbsnp: $dbsnp";
# echo "ID: $ID";

# MuSE
# cmd="MuSE call -f $Reference $TumorBam $NormalBam -O $ID.intermediate"
# echo $cmd; eval $cmd
# 
# cmd="MuSE sump -I $ID.intermediate.MuSE.txt -E -D $dbsnp -O $ID.muse_variants.vcf"
# echo $cmd; eval $cmd &


## MUTECT2 

if [ $gnomad == "null" ]; then gnomad=""; else gnomad="--germline-resource $gnomad";fi
if [ $commonVCF == "null" ]; then commonVCF=""; else commonVCF="-V $commonVCF -L $commonVCF";fi

cmd="gatk Mutect2 -R $Reference -I $TumorBam -I $NormalBam -O $ID.mutect_unfilter.vcf"
cmd="$cmd -L $(echo $Reference|sed 's/fa/interval_list/g') $gnomad -pon pon.vcf.gz"
cmd="$cmd -normal $(samtools view -H $NormalBam|grep ^@RG|sed "s/\t/\n/g"|grep SM|sed 's/SM://g')"
cmd="$cmd --f1r2-tar-gz $ID.f1r2.tar.gz --verbosity ERROR 2>/dev/null"
echo $cmd; eval $cmd

cmd="gatk LearnReadOrientationModel -I $ID.f1r2.tar.gz -O $ID.orientation-model.tar.gz 2>/dev/null"
echo $cmd; eval $cmd

cmd="gatk GetPileupSummaries -I $TumorBam $commonVCF -O ${ID}_T.pileup.tab 2>/dev/null"
echo $cmd; eval $cmd&

cmd="gatk GetPileupSummaries -I $NormalBam $commonVCF -O ${ID}_N.pileup.tab 2>/dev/null"
echo $cmd; eval $cmd&

wait

cmd="gatk CalculateContamination -I ${ID}_T.pileup.tab -matched ${ID}_N.pileup.tab"
cmd="$cmd --tumor-segmentation $ID.segments.tab -O $ID.contamination.tab 2>/dev/null"
echo $cmd; eval $cmd

cmd="gatk FilterMutectCalls -R $Reference -V $ID.mutect_unfilter.vcf -O $ID.mutect_filtered.vcf"
cmd="$cmd --tumor-segmentation $ID.segments.tab --stats $ID.mutect_unfilter.vcf.stats"
cmd="$cmd --contamination-table $ID.contamination.tab 2>/dev/null"
echo $cmd; eval $cmd

cp $ID.mutect_filtered.vcf $ID.final.vcf

# VARSCAN

# cmd="samtools mpileup -f $Reference $NormalBam $TumorBam -q 1 -B > ${ID}.pileup 2>/dev/null"
# echo $cmd; eval $cmd
# 
# cmd="varscan somatic"
# cmd="$cmd ${ID}.pileup $ID.varscan"
# cmd="$cmd --mpileup 1 --min-coverage 8 --min-coverage-normal 8 --min-coverage-tumor 6"
# cmd="$cmd --min-var-freq 0.10 --min-freq-for-hom 0.75 --mpileup 1"
# cmd="$cmd --min-coverage 8 --min-coverage-normal 8 --min-coverage-tumor 6"
# cmd="$cmd --min-var-freq 0.10 --min-freq-for-hom 0.75"
# cmd="$cmd --normal-purity 1.0 --tumor-purity 1.00"
# cmd="$cmd --p-value 0.99 --somatic-p-value 0.05"
# cmd="$cmd --strand-filter 0 --output-vcf"
# echo $cmd; eval $cmd
# 
# cmd="varscan processSomatic $ID.varscan.snp.vcf"
# cmd="$cmd --min-tumor-freq 0.10 --max-normal-freq 0.05 --p-value 0.07"
# echo $cmd; eval $cmd &
# 
# cmd="varscan processSomatic $ID.varscan.indel.vcf"
# cmd="$cmd --min-tumor-freq 0.10 --max-normal-freq 0.05 --p-value 0.07"
# echo $cmd; eval $cmd &
# 
# wait
# 
# cmd="gatk GatherVcfs --INPUT $ID.varscan.snp.vcf --INPUT $ID.varscan.indel.vcf --OUTPUT $ID.varscan.vcf"
# cmd="$cmd --CREATE_INDEX false"
# echo $cmd; eval $cmd

## Merge Mutation Call

# cmd="gatk MergeVcfs --VERBOSITY ERROR --CREATE_INDEX false"
# # cmd="$cmd -I $ID.muse.vcf -I $ID.varscan.vcf -I $ID.mutect_filtered.vcf -O $ID.final.vcf"
# cmd="$cmd $ID.mutect_filtered.vcf -O $ID.final.vcf"
# echo $cmd; eval $cmd


## Strelka2 (need python2)
# cmd="configureStrelkaSomaticWorkflow.py --tumorBam $TumorBam --normalBam $NormalBam"
# cmd="$cmd --referenceFasta $Reference --runDir $ID"
# echo $cmd; eval $cmd
# 
# $ID/runWorkflow.py -m local -j 20
# 
# picard SortVcf 
# 
# cmd="gatk GatherVcfs --INPUT $ID/results/variants/somatic.indels.vcf.gz"
# cmd="$cmd --INPUT $ID/results/variants/somatic.snvs.vcf.gz"
# cmd="$cmd --OUTPUT $ID.Strelka2.vcf --CREATE_INDEX false"
# echo $cmd; eval $cmd


# PINDEL 
# Step 1: Filter Reads

# cmd="Sambamba view $TumorBam --output-filename $TumorBam"
# cmd="$cmd --filter 'not (unmapped or duplicate or secondary_alignment or failed_quality_control or supplementary)'"
# cmd="$cmd --format bam --nthreads 1"
# echo $cmd; eval $cmd


# Step 2: Pindel
# python Codes/PindelSetting.py

# cmd="samtools view -f66 InputData/PBL2_CNV.bin.aligned.sorted.bam | head -n 100000|awk {'print $8'} > $ID.count"
# cmd="cat $ID.count|egrep ^[0-9]{4}$ > $ID."
# printf "$TumorBam\t150\t$ID" > $ID.config_file
# cmd="pindel -f $Reference -i $ID.config_file -o $ID"
# echo $cmd; eval $cmd
# 
# cmd="pindel2vcf -p $ID -r $Reference "

# perl somatic_indelfilter.pl $(somatic.indel.filter.config)

# cmd="picard.jar SortVcf SEQUENCE_DICTIONARY=GRCh38.d1.vd1.dict I=$ID.pindel.somatic.vcf OUTPUT=$ID.pindel.somatic.vcf.gz"
# cmd="$cmd CREATE_INDEX=true"
# echo $cmd; eval $cmd

# cmd="gatk VariantFiltration -V input.vcf.gz -O output.vcf.gz -R GRCh38.d1.vd1.fa"
# cmd="$cmd --disable_auto_index_creation_and_locking_when_reading_rods"
# cmd="$cmd --filterExpression 'vc.isBiallelic() && vc.getGenotype(\"TUMOR\").getAD().1 < 3'"
# cmd="$cmd --filterName TALTDP"
# echo $cmd; eval $cmd

