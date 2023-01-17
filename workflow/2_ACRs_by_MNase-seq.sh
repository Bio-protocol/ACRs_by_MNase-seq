## Identification of accessible chromatin regions with MNase-seq ##
## usage: sh ACRs_by_MNase-seq.sh [name] ##

# 1. Working_path
cd data1/user/path
mkdir MNase_seq
cd MNase_seq
mkdir raw_data

# 2. FastQC_raw
mkdir ./fastqc_result
for i in ./raw_data/*.fq.gz
do
fastqc -o ./fastqc_result $i
done

# 3. Trimmomatic
mkdir clean_data
java -jar trimmomatic-0.36.jar PE -phred33 \
  ./raw_data/*.r1.fq.gz \
  ./raw_data/*.r2.fq.gz \
  ./clean_data/$name.1P.fq \
  ./clean_data/$name.1U.fq \
  ./clean_data/$name.2P.fq \
  ./clean_data/$name.2U.fq \
  ILLUMINACLIP:data1/user/tools/Trimmomatic-0.36/adapters/\
  TruSeq3-PE-2.fa:2:30:10:1:true \
  LEADING:5 TRAILING:5 MINLEN:20 >trimmomatic.$name.log

# 4. FastQC_clean
for i in ./clean_data/*.fq.gz
do
fastqc -o ./fastqc_result $i
done

# 4. Bowtie2
bowtie2-build ref_genome.fa ref.genome
samtools faidx ref_genome.fa
bowtie2 --phred33 -p 10 --reorder -5 6 \
  --no-mixed --no-discordant --no-unal --dovetail \
  -x /index_path/ref.genome \
  -1 ./$name.1P.fq \
  -2 ./$name.2P.fq \
  -S ./$name.sam >bowtie2.$name.log

# 5. Samtools
function samtools_scripts(){
 samtools view -q 20 -bhS $name.sam -o $name.Q20.bam
 samtools sort -n $name.Q20.bam -o $name.Q20.nsort.bam
 samtools fixmate -m $name.Q20.nsort.bam $name.Q20.nsort.ms.bam
 samtools sort $name.Q20.nsort.ms.bam -o $name.Q20.ms.psort.bam
 samtools markdup -r $name.Q20.ms.psort.bam $name.Q20.psort.markdup.bam
 samtools index -c $name.Q20.psort.markdup.bam
}
samtools_scripts $name >samtools.$name.log
rm $name.Q20.bam
rm $name.Q20.nsort.bam
rm $name.Q20.nsort.ms.bam
rm $name.Q20.ms.psort.bam

# 6. "DNS" score
bamCompare -b1 leaf_light_rep1.bam -b2 leaf_heavy_rep1.bam \
  --scaleFactorsMethod None --normalizeUsing CPM --operation subtract \
  --binSize 10 --smoothLength 50 --outFileFormat bedgraph \
  -o leaf_diff_rep1.bdg

# 7. Bayes factor
bamCoverage -b leaf_light_rep1.bam \
  --binSize 10 --smoothLength 50 \
  --normalizeUsing CPM --outFileFormat bedgraph \
  -o leaf_light_rep1.bdg
bamCoverage -b leaf_heavy_rep1.bam \
  --binSize 10 --smoothLength 50 \
  --normalizeUsing CPM --outFileFormat bedgraph \
  -o leaf_heavy_rep1.bdg
bedtools unionbedg \
  -i leaf_light_rep1.bdg leaf_heavy_rep1.bdg \
  -header -names light heavy \
>leaf_rep1.unionbdg
Rscript bayes_factor_caculator.R leaf_rep1.unionbdg

# 8. ACRs
bedtools unionbedg \
  -i leaf_diff_rep1.bdg leaf_rep1.unionbdg.bayes \
  >leaf.diff_rep1_bayes
cat leaf.diff_rep1_bayes | awk '$4>0 && $5>0.5' | \
  cut -f1-3 | bedtools merge -d 200 \
  >leaf.MSF_rep1_bayes_0.5_merge_200.bed
cat leaf.diff_rep1_bayes | awk '$4<0 && $5>0.5' | \
  cut -f1-3 | bedtools merge -d 200 \
  >leaf.MRF_rep1_bayes_0.5_merge_200.bed

# 9. IGV
bamCoverage -b leaf_light_rep1.bam \
--binSize 10 --smoothLength 50 \
  --normalizeUsing CPM --outFileFormat bigwig \
  -o leaf_light_rep1.bw
bamCoverage -b leaf_heavy_rep1.bam \
--binSize 10 --smoothLength 50 \
  --normalizeUsing CPM --outFileFormat bigwig \
  -o leaf_heavy_rep1.bw
bamCompare -b1 leaf_light_rep1.bam -b2 leaf_heavy_rep1.bam \
  --scaleFactorsMethod None --normalizeUsing CPM --operation subtract \
  --binSize 10 --smoothLength 50 --outFileFormat bigwig \
  -o leaf_diff_rep1.bw

