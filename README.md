# BSIAB_Population_in_China
 This repository contains the analysis scripts that were written and used in the article of "**National genomic epidemiology and phylodynamics of Acinetobacter baumannii bloodstream isolates in China: the population superseding and toxification**".
 
 The links below the sub-headings lead to the scripts needed for the corresponding steps. Most of the scripts were developed for running on the Ubuntu 18.04.6 LTS server. You may download and adapt the scripts to suit your own requirements.

## 1. Database and software used in this workflow
### Database
- [COG](https://www.ncbi.nlm.nih.gov/research/cog/)
- [CARD](https://card.mcmaster.ca/)
- [VFDB](http://www.mgc.ac.cn/VFs/)
- [BacMet](http://bacmet.biomedicine.gu.se/)
- [pubMLST](https://pubmlst.org/)
- [isfinder](https://isfinder.biotoul.fr/about.php)
- [ICEberg](https://tool2-mml.sjtu.edu.cn/ICEberg3/index.html)
- [PHASTER](http://phaster.ca/)
- [Islandviewer4](https://www.pathogenomics.sfu.ca/islandviewer)

### Software
- [R](https://www.r-project.org/)
- [Perl](https://www.perl.org/)
- [Python3](https://www.python.org/)
- [PRINSEQ-lite](https://github.com/uwb-linux/prinseq)
- [SPAdes](https://github.com/ablab/spades)
- [QUAST](https://github.com/ablab/quast)
- [FastANI](https://github.com/ParBLiSS/FastANI)
- [Prokka](https://github.com/tseemann/prokka)
- [mlst](https://github.com/tseemann/mlst)
- [Kaptive](https://github.com/klebgenomics/Kaptive)
- [Phytools](https://cran.r-project.org/web/packages/phytools/index.html)
- [RAxML](https://evomics.org/learning/phylogenetics/raxml/)
- [Gubbins](https://github.com/nickjcroucher/gubbins)
- [minimap2](https://github.com/lh3/minimap2)
- [BEAST2](https://www.beast2.org/)
>The following code takes the strain named S1 as an example

## 2. Reads trimming
### PRINSEQ-lite
```
prinseq-lite -verbose -fastq S1_clean_1.fq -ns_max_p 10 -min_qual_mean 25 -out_good S1_1_good -out_bad S1_1_bad

prinseq-lite -verbose -fastq S1_clean_2.fq -ns_max_p 10 -min_qual_mean 25 -out_good S1_2_good -out_bad S1_2_bad
```
## 3. Assembly
### SPAdes
`spades.py -1 S1_1_good.fastq -2 S1_2_good.fastq --phred-offset 33 -o S1 -t 30 -m 1024`

## 4. Genome assembly evaluation
### QUAST
`quast -t 10 -o S1_N50 S1.fasta`

## 5. Amino acid identity (ANI) calculation
### FastANI
`fastANI --ql list.txt --rl strainlist.txt --matrix -t 40  -o ani.txt`

## 6. Genome annotation
### Prokka
`prokka --compliant --kingdom Bacteria --genus Acinetobacter --species baumannii --cpus 20 --outdir outdirname --prefix S1 S1.fasta`

### Compare with database
`diamond blastp --db  database_name.dmnd -ungapped-score 95 --min-score 95 --query S1.faa --out test.tab --query-cover 80 --subject-cover 80 --outfmt 6 --evalue 1e-3 --id 80`

## 7. ST identification and diversity analysis
### mlst
`mlst --scheme abaumannii --threads 40 S1.fasta`

### diversity analysis
Diversity indices of ST s were calculated via a custom script available on GitHub (https://github.com/ahmedmagds/rarefaction-curves). 

You can also find and download the script in `rarefaction-curves` directory.

## 8. Identification of surface polysaccharide synthesis sites
### Kaptive
```
kaptive.py -a S1.fna -k reference_database/Acinetobacter_baumannii_k_locus_primary_reference.gbk -o S1.out

kaptive.py -a S1.fna -k reference_database/Acinetobacter_baumannii_OC_locus_primary_reference.gbk -o S1.out
```
A step-by-step tutorial is available [here](https://bit.ly/kaptive-workshop).

## 9. Core SNP Phylogenetic analysis
```
# call SNPs for multiple isolates from the same reference NC_021726.1, then get clean.full.aln file.

### Gubbins
# detect recombination region
run_gubbins.py -f 50 -p gubbins clean.full.aln

# remove recombination region
snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln
# -c only output columns containing exclusively ACGT

### RAxML
# build core SNP tree
raxmlHPC -f a -x 12345 -p 12345 -# 100 -m GTRGAMMAX -s clean.core.aln -n tre

```

## 10. Ancestral state reconstruction of mobile genetic element (MGE) number
### Phytools, R
```
# SNP.tre and mge_statistics.csv can be found in `MGE_ancestral_state` dictionary
# selected Transposon as example
setwd("path/to/work_dictionary")
library(phytools)

tree <- read.tree("SNP.tre")
#the phylogenetic tree built in 9 above.

mge_info = read.csv(file ='mge_statistics.txt', sep="\t", header=T, row.names = 1,check.names=FALSE)
mge_Tn<-as.matrix(mge)[,1] # select one

# estimate ancestral states and compute variances & 95% confidence intervals for each node:
fit<-fastAnc(tree,mge_Tn,vars=TRUE,CI=TRUE)
#print(fit,printlen=10)
# projection of the reconstruction onto the edges of the tree
obj<-contMap(tree,mge_Tn,plot=FALSE)
plot(obj,legend=0.7*max(nodeHeights(tree)),
     fsize=c(0.2,0.9), lwd=1, outline = F, leg.txt="IS")

# fsize, set font size
# outline, logical value indicated whether or not to outline the plotted color bar with a 1 pt line.
# leg.txt, set title of legend.
# ftype="off", don't show leave's name.

# OR set colors manually
obj<-setMap(obj,c("red", "#fffc00", "green", "purple", "blue", "#d7ff00", "black"))
plot(obj,legend=0.5*max(nodeHeights(tree)),
     fsize=c(0.1,0.9), lwd=1, outline = F, leg.txt="Transposon")

```
