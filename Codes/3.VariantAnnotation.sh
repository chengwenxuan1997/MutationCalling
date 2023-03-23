
CachePath="Reference/homo_sapiens_vep_109_GRCh37"
Reference="Reference/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa"

while getopts c:r: flag
do
    case "${flag}" in
        c) CachePath=${OPTARG};;
        r) Reference=${OPTARG}
    esac
done

Version=$(basename $CachePath|sed 's/.*vep//g'|sed 's/\(_[0-9]*\)_\(.*\)/\1/g'|sed 's/_//g')
Genome=$(basename $CachePath|sed 's/.*vep//g'|sed 's/\(_[0-9]*\)_\(.*\)/\2/g'|sed 's/_//g')

# echo "ID: $ID";
# echo "CachePath: $CachePath";
# echo "Reference: $Reference";
# echo "Genome: $Genome";
# echo "Version: $Version";

# merge all sample
mv $(ls|grep final.vcf) Results
for vcf in $(ls Results|grep final.vcf)
do
{
  TumorID=$(cat Results/$vcf |grep tumor_sample|sed 's/##tumor_sample=//g')
  NormalID=$(cat Results/$vcf |grep normal_sample|sed 's/##normal_sample=//g')
  
  # vep (Variant Effect Predictor)
  cmd="vep -i Results/$vcf -o Results/$(echo $vcf|sed 's/final/anno/g')"
  cmd="$cmd --dir $CachePath --fasta $Reference --assembly $Genome --cache_version $Version"
  cmd="$cmd --cache --offline --fork 4 --vcf --force_overwrite"
  echo $cmd; eval $cmd
  
  # vcf2maf
  cmd="perl Packages/vcf2maf-1.6.21/vcf2maf.pl --inhibit-vep --ref-fasta $Reference"
  cmd="$cmd --input-vcf Results/$(echo $vcf|sed 's/final/anno/g')"
  cmd="$cmd --output-maf Results/$(echo $vcf|sed 's/final.vcf/anno.maf/g')"
  cmd="$cmd --tumor-id $TumorID --normal-id $NormalID"
  echo $cmd; eval $cmd

}
done

cmd=$(ls|grep .final.vcf|sed ':a;N;$!ba; s/\n/ -V /g')
cmd="gatk CombineGVCFs -R $Reference -V $cmd -O Results/Merged.vcf 2>/dev/null"
echo $cmd; eval $cmd

cat Results/*.anno.maf | egrep "^#|^Hugo_Symbol" | head -2 > Results/Merged.maf
cat Results/*.anno.maf | egrep -v "^#|^Hugo_Symbol" >> Results/Merged.maf

