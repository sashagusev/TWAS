---
layout: default
---

## Transcriptome-Wide Association Study (via imputation)

**All code and datasets can be downloaded at [https://data.broadinstitute.org/alkesgroup/TWAS](https://data.broadinstitute.org/alkesgroup/TWAS)**

This page provides the documentation and data for computing association between a molecular phenotype (such as gene expression) and a trait from a large GWAS using summary-level data. We provide precomputed gene expression statistics from three studies to facilitate analysis of GWAS data. Please cite the following manuscript when using this software:

| Integrative approaches for large-scale transcriptome-wide association studies
| Alexander Gusev, Arthur Ko, Huwenbo Shi, Gaurav Bhatia, Wonil Chung, Brenda WJ Penninx, Rick Jansen, Eco JC de Geus, Dorret I Boomsma, Fred A Wright, Patrick F Sullivan, Elina Nikkola, Marcus Alvarez, Mete Civelek, Aldons J Lusis, Terho Lehtimaki, Emma Raitoharju, Mika Kahonen, Ilkka Seppala, Olli Raitakari, Johanna Kuusisto, Markku Laakso, Alkes L Price, Paivi Pajukanta, Bogdan Pasaniuc, 2016 *Nature Genetics*

The IMPG method was developed by Huwenbo Shi as described in:

| Fast and accurate imputation of summary statistics enhances evidence of functional enrichment. 
| Pasaniuc B, Zaitlen N, Shi H, Bhatia G, Gusev A, Pickrell J, Hirschhorn J, Strachan DP, Patterson N, Price AL. 2014 *Bioinformatics*

See also [http://bogdan.bioinformatics.ucla.edu/software/impg/](http://bogdan.bioinformatics.ucla.edu/software/impg/) for details on IMPG

For questions or comments, contact Alexander Gusev [[agusev@hsph.harvard.edu](mailto:agusev@hsph.harvard.edu)] or Bogdan Pasaniuc [[pasaniuc@ucla.edu](mailto:pasaniuc@ucla.edu)].

## Analysis

### Contents

This package contains the binary files and scripts for computing transcriptome-wide association study (TWAS) statistics from pre-computed expression weights, as well as scripts for computing expression weights from raw data.

* The `bin/TWAS.sh` script executes the association pipeline using pre-computed weights.
* The `bin/TWAS_get_weights.sh` script executes the pipeline for computing expression weights from PLINK-format genotype data

The scripts rely on having R installed (and callable by `R`), as well as two pre-compiled binaries.

### Typical use: Pre-computed gene expression weights

This analysis uses the precomputed expression weights and requires only GWAS summary data. A separate package contains the precomputed expression weights for each cohort. The example below 
shows how to run the TWAS with pre-computed weights from the CCDC101 gene and summary data (Z-scores) from the Locke et al. GWAS. To run on your own summary data, simply 
substitute the reference to "WEIGHTS_YFS/CCDC101/CCDC101.LOCKE_BMI.GWAS.zscore" to your own summary data for this gene, in the same format. To run the full TWAS, do so 
for every gene present in the WEIGHTS_*/ directory.

Extract this file (and any other WEIGHTS*bz2 files):

~~~
tar -xjf WEIGHTS_YFS.tar.bz2
~~~

This will create a WEIGHTS_YFS/ directory with subdirectories corresponding to each significantly heritable gene in the YFS cohort. You can then execute a single gene association call as follows:

~~~
bash bin/TWAS.sh {pointer to pre-computed expression weights for this gene} {pointer to Z-score file for this locus} {pointer to output}
~~~

For example, to run this using the GWAS summary statistics from the Locke et al. BMI study on gene CCDC101:

~~~
bash bin/TWAS.sh WEIGHTS_YFS/CCDC101/CCDC101.wgt WEIGHTS_YFS/CCDC101/CCDC101.LOCKE_BMI.GWAS.zscore CCDC101.LOCKE_BMI.TWAS
~~~

This will generate the file `CCDC101.LOCKE_BMI.TWAS.imp`, which contains the imputed Z-score and the precision of the Z-score as second to last and last entries. The final TWAS Z-score is `$5/sqrt($6)`. This file should be identical to the file `WEIGHTS_YFS/CCDC101/CCDC101.LOCKE_BMI.TWAS.imp`

Each directory contains the following files for a given gene:

* CCDC101.wgt.map		: The map for SNPs at the cis locus, with header row indicating format.
* CCDC101.wgt.cor		: The SNP-Expression correlation (i.e. marginal eQTL statistics), one line per SNP.
* CCDC101.wgt.ld		: The SNP-SNP correlation/LD matrix, square.
* CCDC101.LOCKE_BMI.GWAS.zscore : The GWAS summary statistics (from Locke et al. GWAS of BMI) for the cis locus, with header row indicating format.
* CCDC101.LOCKE_BMI.TWAS.imp	: The imputed TWAS Z-score for this locus (i.e. the value we're interested in).

### GWAS summary statistics

The expression weight data includes summary statistics for each locus from the Locke et al. BMI study.
Make sure the summary statistics you use match the strand allele codes as well as the physical positions for the SNPs at each locus. If in doubt, remove any strand ambiguous (A-T, G-C) SNPs.
The summary statistics can be partially overlapping with the expression weight data, but they must have the header and be in sorted physical order.

### Computing expression weights from genotype data
The script for pre-computing expression weights assumes that your data is in a standard PLINK format file (not binary), which contains only the desired SNPs in the cis-locus, and expression set as the phenotype.

The script is executed as follows:

~~~
bash bin/TWAS_get_weights.sh {pointer to prefix for PLINK file} {pointer to output}
~~~

for example, assuming you have a CCDC101.ped and CCDC101.map file (not provided) this would re-generate the expression weights:

~~~
bash bin/TWAS_get_weights.sh CCDC101 WEIGHTS_YFS/CCDC101/CCDC101.wgt
~~~

For maximum power, we recommend computing the genetic value of expression using GEMMA/BSLMM and setting that as the phenotype.

## Reference expression data

The `WEIGHTS_*tar.bz2` contain bzipped tarballs of pre-computed expression weights from three studies. Expression weights were computed from the BLUP/genetic value of expression.

*Warning, the extracted files are large (>1GB). All data is in flat, human readable files.*

### The Netherlands Twin Registry (NTR)
Corresponding to data from Wright et al. 2014 Nat Gen.
1103 significantly heritable genes, computed from 1247 unrelated post-QC samples
expression measured by microarray in blood

### METabolic Sindrome In Men (MET/METSIM)
2146 significantly heritable genes, computed from 563 unrelated post-QC samples
expression measured by RNA-seq in adipose tissue
Note: The exact weights used in the paper are available in old_versions/
weights in this directory include additional samples that were initially excluded due to label-swaps

### Young Finns Study (YFS)
3832 significantly heritable genes, computed from 1264 unrelated post-QC samples
expression measured by microarray in blood

Weights from genes that were significantly heritable in *any* cohort are additionally available in the `ETC/WEIGHTS_NONSIGNIFICANT` directory.
These weights were not used in the paper but are provided for further analyses.

## Omnibus test

We proposed an N degree of freedom test for effect across N tissues after adjusting for predictor correlation (as estimated in the YFS data).

### Reproducing results

To reproduce results from the BMI (Locke et al) TWAS, simply run:

~~~
R --slave < OMNIBUS.R
~~~

The file reads in prediction correlations from `OMNIBUS.CORR.dat` and individual reference TWAS results from `OMNIBUS.ZSCORES.dat`. Both files contain headers describing their contents.

The output looks like the following:

| LOCKE.BMI | AATK | 2 | 0.2120281 |
| LOCKE.BMI | ABCC3 | 2 | 0.8500244
| LOCKE.BMI | ABR | 2 | 0.7403101 
| LOCKE.BMI | ACAP1 | 2 | 0.5919793 
| LOCKE.BMI | ACCS | 3 | 0.1391507 

Where the columns correspond to `<GWAS> <GENE> <# measures/degrees-of-freedom> <P-value>`

### Running your own (with YFS/NTR/MET)

To compute omnibus statistics using your own TWAS results, simply replace 
`OMNIBUS.ZSCORES.dat` with a concatenation of your TWAS Z-scores in the 
same format, where each row is:

~~~
<name of the reference (YFS/NTR/MET)> <GWAS study> <Gene name> <TWAS Z-score>
~~~

The script automatically identifies the duplicate entries for eah GWAS:gene
combination and performs the corresponding N-dof test.

### Running your own (with own weights)

To compute statistics using your own weights and TWAS results, replace
`OMNIBUS.CORR.dat` with the pairwise correlation for your predictors in the
same format, where each row is:

~~~
<expression 1 name> <expression 2 name> <gene name> <signed Pearson correlation>
~~~

These can be computed by predicting from all reference panels into a single cohort
and computing the pearson correlation for every pair of predictors for every gene.
The single cohort should be selected to have  linkage disequilibrium patterns
similar to the TWAS data where the tests will be applied.

The script automatically identifies the number of reference panels and constructs
the corresponding NxN correlation matrix.

### Notes

* We have not evaluated the performance of this test on genes that are only
significantly heritable in a single cohort.

* Because this is a MANOVA, it is much more susceptible than GWAS to highly significant
associations when the inputs were computed incorrectly. For example, individually
weak TWAS Z-scores in opposite directions can become highly significant for predictors
that are expected to be correlated in the same direction. We recommend care when
performing the test (i.e. getting many significant genes is not an indication that 
everything worked correctly) and to double-check that marginal TWAS Z-scores are
at least weakly consistent with an association.