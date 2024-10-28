# Preprocessing Pipeline Script

This script preprocesses a genetic dataset in VCF format, performing multiple steps including sample extraction, file format conversion, filtering by minor allele frequency (MAF) and heterozygosity, kinship calculation, and principal component analysis (PCA).

---

## SLURM Job Parameters

```sh
#!/bin/sh
#SBATCH --ntasks=2                       # Number of tasks (processes) to allocate
#SBATCH --mem-per-cpu=90G                # Memory per CPU
#SBATCH --job-name=preprocess2           # Job name for easier tracking
#SBATCH --time=09:00:00                  # Wall time limit for the job
#SBATCH --partition=batch,jclarke,guest  # Partition to submit the job to
#SBATCH --output=%x.out                  # Standard output file name pattern (%x: job name)
#SBATCH --mail-type=ALL                  # Send email notifications for job start, end, and fail events
#SBATCH --mail-user=jtorres-rodriguez2@unl.edu # Email for job notifications
```

## Loading Required Modules

```sh
ml bcftools/1.10 plink2/2.0a1 gemma/0.98
```

Loads necessary modules:
- **bcftools**: For VCF file manipulation.
- **plink2**: For genetic data manipulation and analysis.
- **gemma**: For kinship and association analysis.

---

## Step-by-Step Process

### 1. Copy the Input VCF File

```sh
cp WiDiv.vcf.gz /scratch/
```

This command copies the compressed VCF file `WiDiv.vcf.gz` to the scratch directory for faster I/O access.

### 2. Extract Samples and Output a Filtered VCF

```sh
bcftools view --samples-file names.693.txt --force-samples -I /scratch/WiDiv.vcf.gz -O z -o WiDiv693.vcf.gz
```

**Description**: Filters the original VCF file by keeping only the samples listed in `names.693.txt`, which includes 693 samples, and outputs a compressed VCF file (`WiDiv693.vcf.gz`).

### 3. Convert VCF to PLINK Format

```sh
plink2 --vcf WiDiv693.vcf.gz --double-id --make-bed --out WiDiv693.vcf
```

**Description**: Converts the VCF file to PLINK binary format (`.bed`, `.bim`, `.fam` files) with default IDs, as PLINK is more efficient for downstream analyses.

### 4. Filter by Minor Allele Frequency (MAF)

```sh
plink2 --bfile WiDiv693.vcf --maf 0.05 --make-bed --out WiDiv693.MAF05.vcf
```

**Description**: Filters out variants with a minor allele frequency (MAF) lower than 0.05. The output file (`WiDiv693.MAF05.vcf`) is in PLINK binary format.

### 5. Hardy-Weinberg Equilibrium Test

```sh
plink2 --bfile WiDiv693.MAF05.vcf --hardy
```

**Description**: Computes Hardy-Weinberg Equilibrium (HWE) statistics for each SNP. The output (`plink2.hardy`) includes HWE statistics for quality control purposes.

### 6. Identify and Exclude High Heterozygosity SNPs

```sh
awk '$8 > 0.05 {print $2}' plink2.hardy > high_het_snps.txt
```

**Description**: Identifies SNPs with heterozygosity greater than 0.05 and outputs these SNP IDs to `high_het_snps.txt`.

### 7. Remove High Heterozygosity SNPs

```sh
plink2 --bfile WiDiv693.MAF05.vcf --exclude high_het_snps.txt --make-bed --out WiDiv693.MAF05.HET0.05.vcf
```

**Description**: Excludes SNPs identified in the previous step (`high_het_snps.txt`) and generates a new PLINK binary file (`WiDiv693.MAF05.HET0.05.vcf`).

### 8. Recalculate Hardy-Weinberg Statistics

```sh
plink2 --bfile WiDiv693.MAF05.HET0.05.vcf --hardy --out plink2.hardy2
```

**Description**: Recalculates HWE statistics after filtering for high heterozygosity SNPs. Outputs are saved in `plink2.hardy2`.

### 9. Calculate Kinship Matrix with GEMMA

```sh
gemma -bfile WiDiv693.MAF05.vcf -gk 1 -o WiDiv693.MAF05.kin.vcf -p pheno_693_values.txt
```

**Description**: Uses GEMMA to calculate the kinship matrix (genetic relationship matrix) for the filtered dataset. The kinship matrix is used for population structure and relatedness correction in association studies.

### 10. Perform Principal Component Analysis (PCA)

```sh
plink2 --bfile WiDiv693.MAF05.vcf --pca 10 --out WiDiv693.MAF05.pca
```

**Description**: Performs PCA to identify population structure in the genetic data, extracting the top 10 principal components.

### 11. Extract Top 5 Principal Components

```sh
awk '{print $3,$4,$5,$6,$7}' WiDiv693.MAF05.pca.eigenvec > WiDiv693.MAF05.5.pca.txt
```

**Description**: Extracts the top 5 principal components from the PCA results file (`WiDiv693.MAF05.pca.eigenvec`) and saves them in `WiDiv693.MAF05.5.pca.txt`.

---

## Output Files

- `WiDiv693.vcf.gz`: Filtered VCF for the selected samples.
- `WiDiv693.MAF05.vcf`: PLINK binary files with MAF filtering.
- `plink2.hardy`: HWE statistics.
- `high_het_snps.txt`: SNPs with heterozygosity > 0.05.
- `WiDiv693.MAF05.HET0.05.vcf`: PLINK binary files with high heterozygosity SNPs removed.
- `plink2.hardy2`: HWE statistics after heterozygosity filtering.
- `WiDiv693.MAF05.kin.vcf`: Kinship matrix.
- `WiDiv693.MAF05.pca`: PCA results.
- `WiDiv693.MAF05.5.pca.txt`: Top 5 principal components.

---

## Notes

- Modify SLURM parameters based on system and project requirements.
- Ensure `names.693.txt` and `pheno_693_values.txt` are correctly formatted and located in the working directory.
