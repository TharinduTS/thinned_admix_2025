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
Then I resturctured files as needed for the compiler by running following command where directories for different runs are. Please disregard 'mkdir: Logs: File exists'
```
mkdir for_compiler
cd for_compiler

for i in {2..5}; do mkdir $i;cd $i;for j in {1..20};do mkdir $j;cp ../../run${j}/*$i* $j;mkdir Logs;mv $j/log${i}.out Logs/K${i}_${j}.log;done;cd ..;done
```
Then I made a directory inside the for_compiler directory for additional needed files
```
mkdir extra_files
```
and downloaded .bed file of my original data from where I converted it to plink to this directory
```
cd extra_files
scp premacht@cedar.computecanada.ca:/scratch/premacht/silurana_papaer_pop_gen_2025/admixture/trop_WGS_outs.bed .
```
created a tab seperated txt file for .ind file in the same directory with individual ID, population and sex (male-1, female-2, unknown-0 as in the following example

Data.ind.txt
```
F_Ghana_WZ_BJE4687_combined__sorted.bam	2	GE
F_IvoryCoast_xen228_combined__sorted.bam	2	IC
F_SierraLeone_AMNH17272_combined__sorted.bam	2	SL
F_SierraLeone_AMNH17274_combined__sorted.bam	2	SL
JBL052_concatscafs_sorted.bam	2	NG
M_Ghana_WY_BJE4362_combined__sorted.bam	1	GE
M_Ghana_ZY_BJE4360_combined__sorted.bam	1	GW
M_SierraLeone_AMNH17271_combined__sorted.bam	1	SL
M_SierraLeone_AMNH17273_combined__sorted.bam	1	SL
XT10_WZ_no_adapt._sorted.bam	2	GE
XT11_WW_trim_no_adapt_scafconcat_sorted.bam	2	GE
XT1_ZY_no_adapt._sorted.bam	1	GE
XT7_WY_no_adapt__sorted.bam	1	GE
all_ROM19161_sorted.bam	2	LB

```
then ran the following bash script saved in for_compiler folder after changing following section in this script accordingly

fn0="/PATH/TO/ADMIXTURE/OUTPUT" # This should be a directory with one folder per K (which itself contains one folder per replicate and one Logs folder with the slurm logs).



## Set ADMIXTURE IndFile (Eigenstrat), ADMIXTURE input .bed file and minimum and maximum K values.
IndFile="/PATH/TO/YOUR/PLINK/DATA/Data.ind" # An Eigenstrat .ind version of the data you converted to plink
bedFile="/PATH/TO/YOUR/PLINK/DATA/Data.pruned.bed" # The input .bed file you ran ADMIXTURE on.
Kmin=2 #The minimum number of Ks you ran
Kmax=5 #The maximum number of Ks you ran

```bash
#!/usr/bin/env bash

## set the path to path with this script on

cd "$(dirname "$0")"

###################################################################
## Compiling ADMIXTURE run CV errors and Q matrices for plotting ##
###################################################################

fn0="." # This should be a directory with one folder per K (which itself contains one folder per replicate and one Logs folder with the slurm logs).
cd ${fn0}
mkdir -p ${fn0}/Plotting


## Set ADMIXTURE IndFile (Eigenstrat), ADMIXTURE input .bed file and minimum and maximum K values.
IndFile="./extra_files/Data.ind.txt" # An Eigenstrat .ind version of the data you converted to plink
bedFile="./extra_files/trop_WGS_outs.bed" # The input .bed file you ran ADMIXTURE on.
Kmin=2 #The minimum number of Ks you ran
Kmax=5 #The maximum number of Ks you ran


## Compile CV Errors
touch CVErrors.txt
for i in $(seq ${Kmax} -1 ${Kmin}); do #seq doesn’t like backwards counting on the clusters, so giving the increment of "-1" fixes the problem.
  (echo $i; grep CV $i/Logs/* | cut -f 4 -d " ") | paste  -d " " - CVErrors.txt > temp_CVErrors
  mv temp_CVErrors CVErrors.txt
done
while read r; do echo ${r% } >>temp_CVErrors; done <CVErrors.txt
mv temp_CVErrors CVErrors.txt
mv CVErrors.txt Plotting/CVErrors.txt


## Compile list of replicates with highest Likelihood
for i in $(seq ${Kmin} ${Kmax}); do grep -H ^Logli $i/Logs/*.log | sort -nrk2 | head -n1; done >Plotting/best_runs.txt


## Compile Q matrices and add ind/pop labels
unset runs
for i in $(seq ${Kmin} ${Kmax}); do
  X=$(grep "K${i}\_" Plotting/best_runs.txt | cut -d "K" -f2 | cut -d ":" -f1 )
  K=$(echo ${X} | cut -f1 -d "_"); Rep=$(echo ${X} | cut -d "_" -f 2)
  runs+="$K/${Rep%.log}/$(basename ${bedFile} .bed).$K.Q "
done
paste -d " " ${runs} >Plotting/temp_data.txt


## Create compiled Q matrix header.
for i in $(seq ${Kmin} ${Kmax}); do
  for x in $(seq 1 ${i}); do
    echo -n "${i}:${x}" >>Plotting/temp_header.txt
    echo -n " " >>Plotting/temp_header.txt
  done
done


## Remove trailing space from header.
echo "" >> Plotting/temp_header.txt
while read r; do echo ${r% } > Plotting/temp_header.txt; done <Plotting/temp_header.txt


## Create individual and Pop list
echo "Ind Pop" >Plotting/temp_pop_labels.txt; awk '{print $1,$3}' $IndFile >>Plotting/temp_pop_labels.txt


## Put together header, data and Pop labels to create compiled data table
cd ${fn0}/Plotting
cat temp_header.txt temp_data.txt >temp_compound.data.txt
paste -d " " temp_pop_labels.txt temp_compound.data.txt  >compound.labelled.QperK.txt 

## Clean up
rm temp_*
```
and ran it
```
bash CompileData.sh
```
Then I saved the R script AdmixturePlotter.R in the for_compiler directory

Made a population order file

poporder.txt
```
SL
LB
IC
GE
GW
NG
```
And ran it changing overriding section as needeed through R without bash
AdmixturePlotter.R
```
#!/usr/bin/env Rscript

## Define functions ----------------------------

## A function that calculates the correlation matrix for a K and it's K-1
correlate_components <- function(k, k_min) {
  start_this_k = 2+sum(1:k-1)-(sum(1:k_min-1))+1
  end_this_k = start_this_k+k-1
  end_prev_k = start_this_k-1
  start_prev_k = end_prev_k-(k-2)
  # print (paste0("This K: ",start_this_k,":",end_this_k))
  # print (paste0("Prev K: ",start_prev_k,":",end_prev_k))
  cor(raw_data[start_prev_k:end_prev_k], raw_data[start_this_k:end_this_k])
}

## A function that returns the column name of the best correlated component in 
##   the K-1 run, for each component in the K run.
fix_colours <- function(k, k_min) {
  ## If the K being processed is the minimum K in the data, keep columns 
  ##  unchanged (since there is nothing to compare to).
  if (k == k_min) {
    return(paste0(k, ":", 1:k))
  }
  
  ## Calculate a correlation matrix for each component of this K to the 
  ##   components of K-1.
  cor_mat <- correlate_components(k, k_min)
  
  ## Find the most correlated component from the last K for each component in 
  ##   this K.
  component_order <- c()
  for (x in 1:k-1) {
    component_order <- append(component_order, which.max(cor_mat[x, ]))
  }
  
  ## If one component in the last K is the most correlated with two components 
  ##   on this K, find the next best correlated component, and if that is 
  ##   unique, assign that as the correct component.
  ## If it is not unique, repeat the process until a unique component is found.
  if (any(duplicated(component_order)) == T && 
      sum(duplicated(component_order)) == 1) {
    duplicate <- which(duplicated(component_order))
    condition <- T
    top_correlates <- c()
    while (condition == T) {
      ## 
      ## until one that isn't already crrelated with another component is found.
      top_correlates <- c(top_correlates, which.max(cor_mat[duplicate, ]))
      cor_mat[duplicate, top_correlates] <- NA
      component_order[duplicate] <- which.max(cor_mat[duplicate, ])
      if (any(duplicated(component_order)) == F) {condition = F}
    }
  } else if (sum(duplicated(component_order)) > 1) {
    stop (paste0("Correlation of components failed. Usually this is caused by high CV errors for some of the components you are trying to plot. 
Please consider limiting your input dataset to K=",k_min," to ",k-1,".

You can use this command to extract the suggested columns from the input file:
    cut -d ' ' -f 1-",sum(seq(k_min,k_max-1))+2,"
    "), call.=FALSE)
  }
  ## If a component hasn't been resolved yet, add it as the newest component.
  missing_component = setdiff(1:k, component_order)
  component_order<-append(component_order, missing_component)
  return(paste0(k,":", component_order))
}

pick_colour <- function(x) {
  return(colours[x])
}
#### MAIN ####

## Load libraries -----------------------------
library(optparse)
library(ggplot2)
library(dplyr, warn.conflicts = F)
library(tidyr)
library(stringr)
library(readr)

## Parse arguments ----------------------------
parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input"), type = 'character', 
                     action = "store", dest = "input", 
                     help = "The input data file. This file should contain all 
                     components per K per indiviual for all K values.")
parser <- add_option(parser, c("-c", "--colourList"), type = 'character',
                     action = "store", dest = "colourList", 
                     help = "A file of desired colours, in R compatible formats.
                     One colour per line.")
parser <- add_option(parser, c("-p", "--popOrder"), type = "character",
                     action = 'store', dest = 'popOrder', 
                     help = "A file containing one population per line in the 
                     desired order.")
parser <- add_option(parser, c("-o", "--outputPlot"), type = "character", 
                     action = 'store', default = "OutputPlot", dest = 'output', 
                     help = "The desired name of the output plot. 
                     [Default: '%default.pdf']")
parser <- add_option(parser, c("-r", "--remove"), type = "logical", 
                     action = 'store_true', default = F, dest = 'remove', 
                     help = "If an order list is provided, should populations not 
                     in the list be removed from the output plot?
                     Not recommended for final figures, but can help in cases 
                     where you are trying to focus on a certain subset of your 
                     populations.")

args <- parse_args(parser)

# ****************THIS SECTION IS COMMENTED OUT BECAUSE i AM OVERRIDING INPUTS*******
## If no input is given, script will exit and provide Usage information.
# if (is.null(args$input) == T) {
#   write("No input file given. Halting execution.", stderr())
#   print_help(parser)
#   quit(status = 1)
# }

#************************************************************************************

## Read cli options into variables.
input <- args$input
colour_file <- args$colourList
## Output name will ignore '.pdf' suffix if provided by user
output <- sub(x = args$output, replacement = "", pattern = ".pdf") 
pop_order <- args$popOrder
if (args$remove && is.null(args$popOrder)){
  write("No population order specified. 'remove' option ignored.", stderr())
}

## Load data --------------------------------

#**************OVERRIDING INPUTS HERE**************
# set working directory to current script directory 
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#set input file
input<-"./Plotting/compound.labelled.QperK.txt"
#set colours
colours <-c('red','green','orange','lightblue','darkblue')
#set pop order
pop_order<-"./extra_files/poporder.txt"

## read data
raw_data <- read_delim(input, " ", col_types = cols())

## Infer min and max K values.
k_min <- as.numeric(str_split_fixed(names(raw_data[,3]),":",2)[1])
k_max <- as.numeric(str_split_fixed(names(raw_data[,ncol(raw_data)]),":",2)[1])

## Sort components of each K according to correlation with components of K-1. 
##   This needs to happen per K, otherwise the correlations will not match 
##   beyond the first pair of Ks.
header <- names(raw_data) ## Take column names from original data

## 'Ind' and 'Pop' should always be at the start of the reformatted data.
refcols = c("Ind", "Pop") 

## For each K in the data, use fix_colours to extract the vector of most correlated
##   column names for each Component in the K. Then sort the components of this 
##   K in the raw data.
for (k in k_min:k_max) { 
  refcols <- c(refcols, fix_colours(k,k_min)) ##    
  raw_data <- raw_data[, c(refcols, setdiff(names(raw_data), refcols))] ## 
}

## Finally, fix the column names so that inference of component numbers is correct
names(raw_data) <- header 

## Flatten data to long format
long_data <- gather(raw_data, temp, value, 3:ncol(raw_data))
## Remove raw_data from memory to reduce memory footprint.
rm(raw_data)
## Split K and Component name to separate columns
long_data <- long_data %>% 
  separate(temp, c("K","Component"), sep = ":") %>%
  mutate(K = as.numeric(K), Component = as.numeric(Component))

## If no colour list is provided, use rainbow() to generate the required number 
##   of colours. Otherwise, read the colour definitions into a vector.
#*******************COMMENTING THIS OUT AS THIS WAS OVERRIDDEN ABOVE********
# if (is.null(colour_file) == T){
#   colours = rainbow(k_max)
# } else {
#   colours = read_delim(colour_file, "\n", col_types = cols(), col_names = F)
#   colours <- colours$X1
# }
#***************************************************************************

## Create colour column based on colour vector.
## Each component in each K run is given the colour of the same index as that 
##   component from the colours list.
long_data <- long_data %>% 
  mutate(clr = purrr::map(Component, pick_colour) %>% 
           unlist)

## Set order of Pops
#getting rid of unintentional \r
long_data$Pop<-gsub('\r','',long_data$Pop)

## If no OrderList is provided, then the populations are sorted alphabetically
if (is.null(pop_order) == F) {
  order <- read.delim(pop_order, header = F, col.names = "Pops")
  long_data$Pop_f <- factor(long_data$Pop, levels = order$Pops)
} else {
  long_data$Pop_f <- long_data$Pop
}

## Early testing dataset subset
# temp_data <- filter(long_data, K == 2)

## Create the named vector (dictionary) of colours needed for scale_fill_manual.
## Each colour is mapped to itself.
col <- as.character(long_data$clr)
names(col) <- as.character(long_data$clr)

if (is.null(args$remove) == F) {
  long_data <- drop_na(long_data, 'Pop_f')
}

## Plot data --------------------------------------------

## Plot the value of each component(y) per individual(x). 'clr' is also the 
##   categorical variable, which is ok since each category will be seen once per K.
ggplot(long_data, aes(x = Ind , y = value, fill = clr)) +
          geom_bar(stat = 'identity', width = 1) +
  ## Colour bars by colour vector.
  scale_fill_manual(values = col) +
  ## X scale changed per Pop. 0 multiplicative change, and +1 additive.
  ## Creates the white bars between groups.
  scale_x_discrete(expand = c(0, 0.55)) +
  ## Set Y axis label
  ylab("") +
  theme_minimal() +
  theme(legend.position = "none", ## No legend
        text = element_text(family = "Helvetica",size = 36),
        # axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        ## Rotate and resize X axis strip text (Pop Names)
        strip.text.x = element_text(angle = 45, hjust = 0.5, size = 22),
        ## Rotate Y axis strip text (K value)
        strip.text.y = element_text(angle = 180),
        panel.spacing.x = unit(0, "lines"),
        ## Set white space between K plots.
        panel.spacing.y = unit(0.005, "lines"), 
        ## Remove axis ticks, and axis text (ancestry proportion and sample names)
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        ## Remove X axis title ("Ind") 
        axis.title.x = element_blank(), 
        ## Remove gridlines from plot
        #panel.grid = element_blank()
        ) +
  ## Creates the plot made so far for each K and each Pop.
  ## The plots per POP are then plotted on top of one another to create each K plot. 
  ## The per K plots are plotted below one another.
  facet_grid(K~ Pop_f,
             scales = "free_x",
             space = "free",
             ## switchlabels of the Y-axis so they is plotted to the left (K=).
             switch = "y")
## Saves the plot as a pdf with specified size.
ggsave(filename = paste0(output,".pdf"), 
        limitsize = F,
        width = 25, height = 10,
        units = "cm")

## Silently remove the Rplots.pdf file, if one was created.
if (file.exists("Rplots.pdf") && output !=  "Rplots") {
    invisible(file.remove("Rplots.pdf"))
}
```

