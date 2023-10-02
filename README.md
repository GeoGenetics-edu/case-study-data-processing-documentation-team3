[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/-7_RZisP)

# Day 1


Loop to process the 5 samples:
```
RED='\033[0;36m'
NC='\033[0m'

for i in ~/course/data/day2/fastq/*.gz
do

#Define the prefixes of all of the fastq files
prefix=$(basename $i .fq.gz)


printf "\n${RED}######################################## \n working on $i \n ########################################${NC} \n"
#SEE INSERT LENGTH = 30


printf "\n${RED}############################## \n remove adapers and reads shorter than 30bp (fastp) \n ##############################${NC}\n"
fastp -i ${i} -o ${prefix}.trimmed.fastq -l 30


printf "\n${RED}############################## \n remove duplicates (vsearch) \n ##############################${NC}\n"
vsearch --fastx_uniques ${prefix}.trimmed.fastq --fastqout ${prefix}.vs.fq --minseqlength 30 --strand both
gzip ${prefix}.vs.fq


printf "\n${RED}############################## \n mapping the reads (bowtie2) \n ##############################${NC}\n"
bowtie2 --threads 5 -k 100 -x ~/course/data/shared/mapping/db/aegenomics.db -U ${prefix}.vs.fq.gz --no-unal | samtools view -bS - > ${prefix}.bam


printf "\n${RED}############################## \n sorting and indexing the alignments (samtools) \n ##############################${NC}\n"
samtools sort ${prefix}.bam -o ${prefix}.sorted.bam
samtools index ${prefix}.sorted.bam


printf "\n${RED}############################## \n investigating ancient damage patterns (mapDamage) \n ##############################${NC}\n"
mapDamage -i ${prefix}.sorted.bam -r ~/course/data/shared/mapping/db/aegenomics.db.fasta --no-stats
done
```

We kept the "standard" parameters

# Day 2
