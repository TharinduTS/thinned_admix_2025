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
