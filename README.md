# thinned_admix_2025

working in Beluga/thinned_admixture_2025

First I copied my vcf here
Then removed missing sites and thinned vcf
```
vcftools --vcf trop_WGS_all_20_samples_all_chrs.vcf --max-missing-count 0 --minQ 30 --thin 5000 --recode --recode-INFO-all --out trop_WGS_all_20_samples_all_chrs_no_missing_thinned_5000.vcf
```
converted to plink
```
#!/bin/sh
#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=30gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load nixpkgs/16.09  intel/2016.4 plink/1.9b_5.2-x86_64
plink --vcf ./trop_WGS_no_cal_mello_niger_all_chrs_thinned_5000.vcf.recode.vcf --make-bed --geno 0.999 --out ./trop_WGS_outs --allow-extra-chr --const-fid
```
--geno 0.999 remove all loci where more than 99.9% of genotypes are missing.

This produces a bunch of support files.

ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0 with following lines
```
FILE=trop_WGS_outs
awk '{$1="0";print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim
```
Then I did multiple runs with the following command
```
#!/bin/sh
#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=8gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben
#SBATCH --array=1-20

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

module load admixture
mkdir run${SLURM_ARRAY_TASK_ID}
cd run${SLURM_ARRAY_TASK_ID}
cp ../trop_WGS_outs* .

for i in {2..5}; do  admixture  -C 0.0000001 -s $((1 + $RANDOM % 1000)) --cv trop_WGS_outs.bed $i> log${i}.out; done
```
Then I downloaded only the necessery directories 
```
rsync -axvH --no-g --no-p --include='run*/***' --exclude='*' premacht@cedar.computecanada.ca:/scratch/premacht/silurana_papaer_pop_gen_2025/admixture/ .
```
Then I copied two following files in each directory, ran R script and changed colours accordingly
trop_WGS_outs.list (list of samples made in the same order as initial file
```
F_Ghana_WZ_BJE4687_combined__sorted.bam Ghana
F_IvoryCoast_xen228_combined__sorted.bam IvoryCoast
F_SierraLeone_AMNH17272_combined__sorted.bam SierraLeone
F_SierraLeone_AMNH17274_combined__sorted.bam SierraLeone
JBL052_concatscafs_sorted.bam Ref
M_Ghana_WY_BJE4362_combined__sorted.bam Ghana
M_Ghana_ZY_BJE4360_combined__sorted.bam Ghana
M_SierraLeone_AMNH17271_combined__sorted.bam SierraLeone
M_SierraLeone_AMNH17273_combined__sorted.bam SierraLeone
XT10_WZ_no_adapt._sorted.bam Lab
XT11_WW_trim_no_adapt_scafconcat_sorted.bam Lab
XT1_ZY_no_adapt._sorted.bam Lab
XT7_WY_no_adapt__sorted.bam Lab
all_ROM19161_sorted.bam Liberia
```
plot_admixture_with_admix_2024.R
```
#!/usr/bin/Rscript

# Usage: plotADMIXTURE.r -p <prefix> -i <info file, 2-column file with ind name and population/species name> 
#                        -k <max K value> -l <comma-separated list of populations/species in the order to be plotted>
# This R script makes barplots for K=2 and all other K values until max K (specified with -k). It labels the individuals 
# and splits them into populations or species according to the individual and population/species names in the 2-column file specified with -i.
# The order of populations/species follows the list of populations/species given with -l.
# Usage example: plotADMIXTURE.r -p fileXY -i file.ind.pop.txt -k 4 -pop pop1,pop2,pop3
# In this example, the script would use the files fileXY.2.Q, fileXY.3.Q, fileXY.4.Q to make barplots for the three populations.
# file.ind.pop.txt should contain one line for each individual in the same order as in the admixture files e.g.
# ind1 pop1
# ind2 pop1
# ind3 pop2
# ind4 pop3

# Author: Joana Meier, September 2019
# edited by Tharindu, Jul 2024

library("rstudioapi")
library(fuzzyjoin)
library(dplyr)

#set the directory which your R script is in tas the working directory
setwd(dirname(getActiveDocumentContext()$path))

# Read in the arguments
library("optparse")
option_list = list(
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="prefix name (with path if not in the current directory)", metavar="character"),
  make_option(c("-i", "--infofile"), type="character", default=NULL, 
              help="info text file containing for each individual the population/species information", metavar="character"),
  make_option(c("-k", "--maxK"), type="integer", default=NULL, 
              help="maximum K value", metavar="integer"),
  make_option(c("-m", "--minK"), type="integer", default=2, 
              help="minimum K value", metavar="integer"),
  make_option(c("-l", "--populations"), type="character", default=NULL, 
              help="comma-separated list of populations/species in the order to be plotted", metavar="character"),
  make_option(c("-o", "--outPrefix"), type="character", default="default", 
              help="output prefix (default: name provided with prefix)", metavar="character")
) 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#******* CHANGE THESE OPTIONS HERE *********

opt$prefix='trop_WGS_outs'
opt$infofile='trop_WGS_outs.list'
opt$maxK=5

# Change pop order here
opt$populations='SierraLeone,Liberia,IvoryCoast,Ghana,Lab,Ref'

#*******************************************

# Check that all required arguments are provided
if (is.null(opt$prefix)){
  print_help(opt_parser)
  stop("Please provide the prefix", call.=FALSE)
}else if (is.null(opt$infofile)){
  print_help(opt_parser)
  stop("Please provide the info file", call.=FALSE)
}else if (is.null(opt$maxK)){
  print_help(opt_parser)
  stop("Please provide the maximum K value to plot", call.=FALSE)
}else if (is.null(opt$populations)){
  print_help(opt_parser)
  stop("Please provide a comma-separated list of populations/species", call.=FALSE)
}

# If no output prefix is given, use the input prefix
if(opt$outPrefix=="default") opt$outPrefix=opt$prefix

# Assign the first argument to prefix
prefix=opt$prefix

# Get individual names in the correct order
labels<-read.table(opt$infofile)

# Name the columns
names(labels)<-c("ind","pop")


# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
labels$n<-factor(labels$pop,levels=unlist(strsplit(opt$populations,",")))
levels(labels$n)<-c(1:length(levels(labels$n)))
labels$n<-as.integer(as.character(labels$n))

# read in the different admixture output files
minK=opt$minK
maxK=opt$maxK
tbl<-lapply(minK:maxK, function(x) read.table(paste0(prefix,".",x,".Q")))


#assign colors to populations at once

col1<-"red"
col2<-"green"
col3<-"orange"
col4<-"skyblue"
col5<-"darkblue"
col6<-"yellow"
col7<-"pink"
col8<-"purple"
col9<-"grey"
col10<-"black"
col11<-"forestgreen"
col12<-"brown"

#create color palettes for each k value
col_palette_k2<-c(col1,col2)
col_palette_k3<-c(col3,col2,col1)
col_palette_k4<-c(col2,col1,col4,col3)
col_palette_k5<-c(col4,col5,col3,col1,col2)


col_palette_k6<-c(col1,col2,col6,col5,col4,col3)
col_palette_k7<-c(col1,col7,col5,col6,col4,col3,col2)
col_palette_k8<-c(col5,col8,col1,col7,col6,col3,col4,col2)
col_palette_k9<-c(col8,col6,col7,col2,col5,col4,col3,col1,col9)
col_palette_k10<-c(col8,col7,col4,col10,col2,col9,col6,col1,col3,col5)
col_palette_k11<-c(col6,col10,col2,col4,col1,col7,col9,col8,col3,col11,col5)
col_palette_k12<-c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12)


# Prepare spaces to separate the populations/species

#****** CHANGE THESE*******

# space_between_individuals

between_inds<-0.05
between_pops<-0.2

#*************************

rep<-as.vector(table(labels$n))
spaces<-0
for(i in 1:length(rep)){spaces=c(spaces,rep(between_inds,rep[i]-1),between_pops)}
spaces<-spaces[-length(spaces)]

# Plot the cluster assignments as a single bar for each individual for each K as a separate row
tiff(file=paste0(opt$outPrefix,".tiff"),width = 800, height = 350,res=200)
par(mfrow=c(maxK-1,1),mar=c(0,1,0,0),oma=c(2,1,1,1),mgp=c(0,0.2,0),xaxs="i",cex.lab=1.2,cex.axis=0.8)
# Plot minK
bp<-barplot(t(as.matrix(tbl[[1]][order(labels$n),])), col=get(paste("col_palette_k",minK,sep="")), border=NA,ylab=minK,yaxt="n",space=spaces)

# add or remove sample names by uncommenting following line. IF YOU COMMENT FOLLOWING, CHANGE 3RD VALUE IN 'ama=' ABOVE TO REMOVE EXTRA SPACE ON TOP

#axis(3,at=bp,labels=labels$ind[order(labels$n)],las=2,tick=F,cex=0.6)

#place this after 'at=' below to change label places manually

#change label locations here
x_ticks<-c(2,4.5,6.4,8.6,12,15)

# place this after 'at=' below to change label places automatically
#c(which(spaces==between_pops),bp[length(bp)])-diff(c(1,which(spaces==between_pops)

# Plot higher K values
if(maxK>minK)lapply(2:(maxK-1), function(x) barplot(t(as.matrix(tbl[[x]][order(labels$n),])), col=get(paste("col_palette_k",x+1,sep="")),xaxt="n", border=NA,ylab=x+1,yaxt="n",space=spaces))
axis(1,at=x_ticks,
     labels=unlist(strsplit(opt$populations,",")))
dev.off()

```
