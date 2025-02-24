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
rsync -axvH --no-g --no-p --include='run*/***' --exclude='*' premacht@cedar.computecanada.ca:/scratch/premacht/silurana_papaer_pop_gen_2025/admixture/ .```
