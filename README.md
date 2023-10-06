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
## Taxonomic profiling and DNA damage estimates (metaDMG) for each individual taxonomic node
We used metaDMG (https://github.com/metaDMG-dev/metaDMG-core) to taxonomically classify each individual read to the loest taxonomic node possible using a similarity identity between 95-100% to reference.
```
conda activate metaDMG

metaDMG config *.sorted.bam --names ~/course/data/shared/mapping/taxonomy/names.dmp --nodes ~/course/data/shared/mapping/taxonomy/nodes.dmp --acc2tax ~/course/data/shared/mapping/taxonomy/acc2taxid.map.gz -m /usr/local/bin/metaDMG-cpp

metaDMG compute config.yaml

metaDMG compute config.yaml 
```


### metaDMG [results](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team3/blob/main/config.yaml)  for eukaryotic taxa at the genus level:

<p align="center">
  <img src="https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team3/assets/48062644/8ce1f1a7-667f-42aa-8e27-c4b97496cb88">
</p>  

We found a shift in the damage distribution between the different samples, which could be an indicator of different sampling depths and, thereby, different ages of the samples. This pattern is clear for the eukaryotic taxa, whereas for bacteria and archea, no pattern can be seen (and the majority of the latter have poor damage patterns). We speculate that this could be due to the differences in the living conditions for those taxa - where some bacteria could live in the deeper soil levels and most of the eukaryotic deposits come from stationary sources. 
The red sample (plot) having been samples at the deepest level, followed by blue, green, purple and grey. 

# Day 3
## Analysis of output from metaDMG


###Load R packages:

```
library(tidyverse) 
library(reshape2)
library(vegan)
library(rioja)
library(ggplot2)
library(dplyr)
library(gghighlight)
library(ggpubr)
```


###import files:

```
df <- read_csv("metaDMGresults.csv")
head(df)
metaDATA <- read.delim("metadata.tsv")
head(metaDATA)
#change colnames
colnames(metaDATA)[colnames(metaDATA) == "sample_name"] <- "sample"
colnames(metaDATA)[colnames(metaDATA) == "years_bp"] <- "YearsBP"
#merge files
dt <- merge(df, metaDATA, by = "sample")
```
### set filter parameters and filter out only Viridiplantae

```
DamMin2 = 0.00
MapSig2 = 0
MinRead2 = 100
MinLength = 35
dt2 <- dt %>% filter(MAP_damage > DamMin2, N_reads >= MinRead2, mean_L > MinLength, MAP_significance  > MapSig2,  grepl("Viridiplantae",tax_path), grepl("\\bgenus\\b", tax_rank), grepl("", sample))
```


### some plots

```
# reorder samples

dt2$sample <- as.character(dt2$sample)
dt2$sample <- factor(dt2$sample ,  levels=c("PRI-TJPGK-CATN-96-98",
 "PRI-TJPGK-CATN-112-114", "PRI-TJPGK-CATN-160-162","PRI-TJPGK-CATN-224-226",
 "PRI-TJPGK-CATN-288-290"))


ggplot() +
  geom_jitter(data = dt2, aes(x=as.numeric(YearsBP), y=MAP_damage, size = N_reads,color=sample), alpha =0.5) +
  gghighlight(N_reads > 500) +
  xlab("Years BP") +
  ylab("DNA damage") +
  labs(title = "Values for taxa with >500 reads", size = "Number of reads")

```

[plot1_color.pdf](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team3/files/12821317/plot1_color.pdf)


```
GC VS DAMAGE
ggplot(dt2, aes(as.numeric(mean_GC),MAP_damage)) + 
geom_point(aes(color=sample))+ 
geom_smooth(method = lm)+
stat_cor()
```
[plotGC_VS_DAMAGE_color.pdf](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team3/files/12821360/plotGC_VS_DAMAGE_color.pdf)
```
significance VS DAMAGE
ggplot(dt2, aes(as.numeric(MAP_significance),MAP_damage)) + 
geom_point(aes(color=sample))+ 
geom_smooth(method = lm)+
stat_cor()+  gghighlight(N_reads > 1000)+
ggtitle("highlight(N_reads > 1000)")
```

[plot3_color_modi.pdf](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team3/files/12821553/plot3_color_modi.pdf)
```
my_comparisons <- list( c("PRI-TJPGK-CATN-96-98", "PRI-TJPGK-CATN-112-114"), 
c("PRI-TJPGK-CATN-96-98", "PRI-TJPGK-CATN-160-162"), 
c("PRI-TJPGK-CATN-96-98", "PRI-TJPGK-CATN-224-226"), 
c("PRI-TJPGK-CATN-96-98", "PRI-TJPGK-CATN-288-290"), 
c("PRI-TJPGK-CATN-112-114", "PRI-TJPGK-CATN-160-162"),
c("PRI-TJPGK-CATN-112-114", "PRI-TJPGK-CATN-224-226") ,
c("PRI-TJPGK-CATN-112-114", "PRI-TJPGK-CATN-288-290") ,
c("PRI-TJPGK-CATN-160-162", "PRI-TJPGK-CATN-224-226") ,
c("PRI-TJPGK-CATN-160-162", "PRI-TJPGK-CATN-288-290"),
c("PRI-TJPGK-CATN-224-226", "PRI-TJPGK-CATN-288-290") )
# 
ggplot(dt2, aes(sample,as.numeric(MAP_damage))) +
geom_boxplot(aes(color=sample),alpha=0.1)+ 
stat_compare_means(comparisons = my_comparisons)+
  gghighlight(N_reads > 1000)+
  ggtitle("highlight(N_reads > 1000)")
```
[boxplot_damage.pdf](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team3/files/12821563/boxplot_damage.pdf)



# Day 4
## Analysis of microbial aDNA

[plot (2).pdf](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team3/files/12821573/plot.2.pdf)

## The difference between reads and genome abundances ^^
[plot2.pdf](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team3/files/12821572/plot2.pdf)

![Screenshot 2023-10-05 192312](https://github.com/GeoGenetics-edu/case-study-data-processing-documentation-team3/assets/48062644/bc9be8e2-00db-44af-be22-01167f46e737)


