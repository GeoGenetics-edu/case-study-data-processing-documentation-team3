[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/-7_RZisP)

# Day 1


Loop to process the 5 samples:
```
conda activate day1

RED='\033[0;36m'
NC='\033[0m'

#SPECIFY THE CORRECT PATH HERE:
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
samtools sort -n ${prefix}.bam -@ 5 -o ${prefix}.sorted.bam
samtools index ${prefix}.sorted.bam


printf "\n${RED}############################## \n investigating ancient damage patterns (metaDMG) - CANCELED \n ##############################${NC}\n"
# mapDamage -i ${prefix}.sorted.bam -r ~/course/data/shared/mapping/db/aegenomics.db.fasta --no-stats
# metaDMG-cpp lca -bam ${prefix}.sorted.bam -names ~/course/data/shared/mapping/taxonomy/names.dmp -nodes ~/course/data/shared/mapping/taxonomy/nodes.dmp -acc2tax ~/course/data/shared/mapping/taxonomy/acc2taxid.map.gz -weighttype 1 -fix_ncbi 0 -out ${prefix}

done
```

We kept the same parameters as in the example commands. 

# Day 2

```
conda activate metaDMG

metaDMG config *.sorted.bam --names ~/course/data/shared/mapping/taxonomy/names.dmp --nodes ~/course/data/shared/mapping/taxonomy/nodes.dmp --acc2tax ~/course/data/shared/mapping/taxonomy/acc2taxid.map.gz -m /usr/local/bin/metaDMG-cpp

metaDMG compute config.yaml

metaDMG compute config.yaml 
```

![newplot (4)](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team3/assets/48062644/8ce1f1a7-667f-42aa-8e27-c4b97496cb88)

