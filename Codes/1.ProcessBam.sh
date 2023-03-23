# Genome Alignment

bam="InputData/KAP2T_CNV.bin.aligned.sorted.bam"
Reference="Reference/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa"
thread=10
id="test"
dbsnp="Reference/GRCh37.p13.dbSNP151.All.vcf.gz"
realignment="F"
class="Tumor"

while getopts b:r:t:i:d:a:c: flag
do
    case "${flag}" in
        b) bam=${OPTARG};;
        r) Reference=${OPTARG};;
        t) thread=${OPTARG};;
        i) id=${OPTARG};;
        d) dbsnp=${OPTARG};;
        a) realignment=${OPTARG};;
        c) class=${OPTARG}
    esac
done


# echo "bam: $bam";
# echo "Reference: $Reference";
# echo "thread: $thread";
# echo "readlength $readlength";
# echo "id $id";
# echo "dbsnp $dbsnp";
# echo "realignment $realignment";

# STEP 1: clean the bam file

if [ $realignment == "T" ];then
  cmd="samtools view -@ $thread -h $bam | sed 's/\([\.\/][12]\)\(\t[0-9]*\)/\2/g' | samtools view -@ $thread -Sb -o $id.bam"
  echo $cmd; eval $cmd
  
  cmd="samtools sort -n -o $id.sorted.bam $id.bam"
  echo $cmd; eval $cmd
  
  cmd="samtools fastq -@ $thread $id.sorted.bam -1 $id.paired1.fq.gz -2 $id.paired2.fq.gz -0 /dev/null -s /dev/null -n"
  # bedtools bamtofastq -i tmp2.bam -fq paired1.fq -fq2 paired2.fq
  echo $cmd; eval $cmd
  
  cmd="bwa mem -t $thread -T 0 -R \"@RG\tID:$id\tSM:$id\tLB:WES\tPL:MD\" $Reference $id.paired1.fq.gz $id.paired2.fq.gz > $id.realign.sam 2>/dev/null"
  echo $cmd; eval $cmd
  
  cmd="samtools view -Shb $id.realign.sam -o $id.realign.bam"
  echo $cmd; eval $cmd
  
  # bwa aln -t 8 <Reference> <fastq_1.fq.gz> > <sai_1.sai> &&
  # bwa aln -t 8 <Reference> <fastq_2.fq.gz> > <sai_2.sai> &&
  # bwa sampe -r <read_group> <Reference> <sai_1.sai> <sai_2.sai> <fastq_1.fq.gz> <fastq_2.fq.gz> | 
  # samtools view -Shb -o <output.bam> -
  # # If the quality scores are encoded as Illumina 1.3 or 1.5, use BWA aln with the "-l" flag.
  
  cmd="samtools sort --thread $thread -o $id.realign.sorted.bam $id.realign.bam"
  echo $cmd; eval $cmd
  
  cmd="picard MarkDuplicates INPUT=$id.realign.sorted.bam OUTPUT=$id.clean.bam METRICS_FILE=$id.marked_dup_metrics.txt"
  cmd="$cmd CREATE_INDEX=true VALIDATION_STRINGENCY=STRICT 2>/dev/null"
  echo $cmd; eval $cmd
  
  ls|grep $id|grep -v clean|xargs rm
fi

if [ $realignment == "F" ];then
  # remove unpaired reads
  # cmd="Sambamba view $bam --output-filename $id.filter.bam"
  # cmd="$cmd --filter 'not (unmapped or duplicate or secondary_alignment or failed_quality_control or supplementary)'"
  # cmd="$cmd --format bam --nthreads 1"
  # echo $cmd; eval $cmd
  
  # remove unnecessary contig 
  # cat $Reference.fai |awk '{print $1}'|sort > $id.r1.txt
  # samtools view -H $bam|grep SN|awk '{print $2}'|sed 's/SN://g'|sort > $id.r2.txt
  # sort $id.r1.txt $id.r1.txt $id.r2.txt |uniq -u|sed 's/^/grep -v /g' | tr "\n" "\|" > $id.r3.txt
  # cmd="samtools view -H $bam| $(cat $id.r3.txt) cat > $id.clean.sam"
  # echo $cmd;eval $cmd
  # cmd="comm -1 -2 $id.r1.txt $id.r2.txt| xargs samtools view --thread $thread $bam >> $id.clean.sam"
  # echo $cmd; eval $cmd
  # samtools view --thread $thread -bh $id.clean.sam > $id.clean.bam
  
  # cp $bam $id.clean.bam
  # remove chr prefix if it exists
  samtools view -H $bam | sed 's/SN:chr\([0-9XY]\)/SN:\1/g; s/SN:chrM/SN:MT/g'| samtools reheader - $bam > $id.clean.bam
  
  cmd="samtools index -@ 10 $id.clean.bam"
  echo $cmd; eval $cmd
fi

# STEP 6: BaseRecalibrator; DBSNP V.144 
cmd="gatk BaseRecalibrator -I $id.clean.bam -O $id.bqsr.grp"
cmd="$cmd -R $Reference -known-sites $dbsnp  --verbosity ERROR 2>/dev/null"
echo $cmd; eval $cmd
  
# STEP 7: PrintReads (ApplyBQSR in gatk4)
cmd="gatk ApplyBQSR -R $Reference -I $id.clean.bam --bqsr-recal-file $id.bqsr.grp -O $id.clean.calib.bam"
cmd="$cmd --verbosity ERROR 2>/dev/null"
echo $cmd; eval $cmd

# Step8: 
if [ $class == "Normal" ];then
  cmd="gatk Mutect2 -R $Reference -I $id.clean.calib.bam --max-mnp-distance 0 -O $id.normal.vcf.gz 2>/dev/null"
  echo $cmd; eval $cmd
fi


# # STEP 1: RealignerTargetCreator 
# java -jar gatk  RealignerTargetCreator \
# -R <Reference>
# -known <known_indels.vcf>
# -I <input.bam>
# -o <realign_target.intervals>
#   
# # STEP 2: INDELREALIGNER 
# java -jar GenomeAnalysisTK.jar -T IndelRealigner \
# -R <Reference> \
# -known <known_indels.vcf> \
# -targetIntervals <realign_target.intervals> \
# --noOriginalAlignmentTags \
# -I <input.bam> \
# -nWayOut <output.map>