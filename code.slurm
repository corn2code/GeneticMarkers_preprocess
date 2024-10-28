#!/bin/sh
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=90G
#SBATCH --job-name=preprocess2
#SBATCH --time=09:00:00
#SBATCH --partition=batch,jclarke,guest
#SBATCH --output=%x.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xx

ml bcftools/1.10 plink2/2.0a1 gemma/0.98

#cp WiDiv.vcf.gz /scratch/

#bcftools view --samples-file names.693.txt --force-samples -I /scratch/WiDiv.vcf.gz -O z -o WiDiv693.vcf.gz

#plink2 --vcf WiDiv693.vcf.gz --double-id --make-bed --out WiDiv693.vcf

#plink2 --bfile WiDiv693.vcf --maf 0.05 --make-bed --out WiDiv693.MAF05.vcf

#gemma -bfile WiDiv693.MAF05.vcf -gk 1 -o WiDiv693.MAF05.kin.vcf -p pheno_693_values.txt

#plink2 --bfile WiDiv693.MAF05.vcf --hardy # this will generate the plink2.hardy file

awk '$8 > 0.05 {print $2}' plink2.hardy > high_het_snps.txt

plink2 --bfile WiDiv693.MAF05.vcf --exclude high_het_snps.txt --make-bed --out WiDiv693.MAF05.HET0.05.vcf

plink2 --bfile WiDiv693.MAF05.HET0.05.vcf --hardy --out plink2.hardy2


#plink2 --bfile WiDiv693.MAF05.vcf --pca 10 --out WiDiv693.MAF05.pca

#awk '{print $3,$4,$5,$6,$7}' WiDiv693.MAF05.pca.eigenvec > WiDiv693.MAF05.5.pca.txt

