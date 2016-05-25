---
layout: default
---

## Transcriptome-Wide Association Study (via imputation)

This page provides the documentation and data for estimating association statistics between gene expression and phenotypes from summary-level GWAS data. We provide precomputed gene expression weights from three studies to facilitate analysis. Please refer to the following manuscripts for details:

| Gusev et al. "Integrative approaches for large-scale transcriptome-wide association studies"  2016 *Nature Genetics*

| Pasaniuc et al. "Fast and accurate imputation of summary statistics enhances evidence of functional enrichment" 2014 *Bioinformatics*

See also [http://bogdan.bioinformatics.ucla.edu/software/impg/](http://bogdan.bioinformatics.ucla.edu/software/impg/) for details on IMPG

For questions or comments, contact Alexander Gusev [[agusev@hsph.harvard.edu](mailto:agusev@hsph.harvard.edu)] or Bogdan Pasaniuc [[pasaniuc@ucla.edu](mailto:pasaniuc@ucla.edu)].

**Due to large size, all code and datasets are hosted at [https://data.broadinstitute.org/alkesgroup/TWAS](https://data.broadinstitute.org/alkesgroup/TWAS)**

## Outline

1. [Typical TWAS analysis](#typical-twas-analysis)
2. [Reference expression data](#reference-expression-data)
3. [Computing your own expression weights](#computing-your-own-expression-weights)
4. [Testing multiple tissues jointly](#testing-multiple-tissues-jointly)
5. [FAQ](#faq)
6. [Change Log](#change-log)

---

## Typical TWAS analysis

The `TWAS.*.tar.bz2` package contains the binary files and scripts for computing TWAS statistics from pre-computed expression weights, as well as scripts for computing expression weights from your own expression + SNP data.

* The `bin/TWAS.sh` script executes the association pipeline using pre-computed weights.
* The `bin/TWAS_get_weights.sh` script executes the pipeline for computing expression weights from PLINK-format genotype data

The scripts rely on having R installed (and callable by `R`), as well as two pre-compiled binaries.

### Pre-computed gene expression weights

This analysis uses the precomputed expression weights and requires only GWAS summary data (Z-scores or effect sizes and standard errors). A separate package contains the precomputed expression weights for each cohort. The example below 
shows how to run the TWAS with pre-computed weights from the CCDC101 gene and summary data (Z-scores) from the Locke et al. GWAS. To run the full TWAS, repeat this procedure for every gene present in the `WEIGHTS_*/` directory.

Extract this file (and any other WEIGHTS*bz2 files):

~~~
tar -xjf WEIGHTS_YFS.tar.bz2
~~~

This will create a `WEIGHTS_YFS/` directory with subdirectories corresponding to each significantly heritable gene in the YFS cohort. You can then execute a single gene association call as follows:

~~~
bash bin/TWAS.sh \
{pointer to pre-computed expression weights for this gene} \
{pointer to Z-score file for this locus} \
{pointer to output}
~~~

For example, to run this using the GWAS summary statistics from the Locke et al. BMI study on gene CCDC101:

~~~
bash bin/TWAS.sh \
WEIGHTS_YFS/CCDC101/CCDC101.wgt \
WEIGHTS_YFS/CCDC101/CCDC101.LOCKE_BMI.GWAS.zscore \
CCDC101.LOCKE_BMI.TWAS
~~~

This will generate the file `CCDC101.LOCKE_BMI.TWAS.imp`, which contains the imputed Z-score and the precision of the Z-score as second to last and last entries. The final TWAS Z-score is computed as `$5/sqrt($6)` from this file. The output should be identical to the file `WEIGHTS_YFS/CCDC101/CCDC101.LOCKE_BMI.TWAS.imp`

Each directory contains the following files for a given gene:

| *.wgt.map	| The map for SNPs at the cis locus, with header row indicating format.
| *.wgt.cor	| The SNP-Expression correlation (i.e. marginal eQTL statistics), one line per SNP.
| *.wgt.ld	| The SNP-SNP correlation/LD matrix, square.
| *.LOCKE_BMI.GWAS.zscore | The GWAS summary statistics (from Locke et al. GWAS of BMI) for the cis locus, with header row indicating format.
| *.LOCKE_BMI.TWAS.imp	| The imputed TWAS Z-score for this locus (i.e. the value we're interested in).

### GWAS summary statistics

The expression weight data includes summary statistics for each locus from the Locke et al. BMI study.
Make sure the summary statistics you use match the strand allele codes as well as the physical positions for the SNPs at each locus. If in doubt, remove any strand ambiguous (A-T, G-C) SNPs.
The summary statistics can be partially overlapping with the expression weight data, but they must have the header and be in sorted physical order.

---

## Reference expression data

**Due to large size, weights can be downloaded from [https://data.broadinstitute.org/alkesgroup/TWAS](https://data.broadinstitute.org/alkesgroup/TWAS)**

The `WEIGHTS_*tar.bz2` contain bzipped tarballs of pre-computed expression weights from three studies. Expression weights were computed from the BLUP/genetic value of expression.

| **The Netherlands Twin Registry (NTR) based on expression from [Wright et al. 2014 Nat Genet]** |
| 1,103 significantly heritable genes, computed from 1,247 unrelated post-QC samples expression measured by microarray in blood |

| **METabolic Sindrome In Men (MET/METSIM)** |
| 2,146 significantly heritable genes, computed from 563 unrelated post-QC samples with expression measured by RNA-seq in adipose tissue. *Note: The exact weights used in the paper are available in the `old_versions/` directory; weights in this directory include additional samples that were initially excluded due to label-swaps.* |

| **Young Finns Study (YFS)** |
| 3,832 significantly heritable genes, computed from 1,264 unrelated post-QC samples expression measured by microarray in blood. |

Weights from genes that were significantly heritable in *any* cohort are additionally available in the `ETC/WEIGHTS_NONSIGNIFICANT` directory.
These weights were not used in the paper but are provided for further analyses.

*Warning, the extracted files are large (>1GB). All data is in flat, human readable files.*

---

## Computing your own expression weights
The script for pre-computing expression weights assumes that your data is in a standard PLINK format file (ped/map, not binary), which contains only the desired SNPs in the cis-locus, and expression set as the phenotype.

The script is executed as follows:

~~~
bash bin/TWAS_get_weights.sh {pointer to prefix for PLINK file} {pointer to output}
~~~

for example, assuming you have a CCDC101.ped and CCDC101.map file (not provided) the following command would re-generate the expression weights:

~~~
bash bin/TWAS_get_weights.sh CCDC101 WEIGHTS_YFS/CCDC101/CCDC101.wgt
~~~

For maximum power, we recommend computing the genetic value of expression using BLUP/BSLMM and setting that as the phenotype.

---

## Testing multiple tissues jointly

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

### NB

* We have not evaluated the performance of this test on genes that are only
significantly heritable in a single cohort.

* Because this is a MANOVA, it is much more susceptible than GWAS to highly significant
associations when the inputs were computed incorrectly. For example, individually
weak TWAS Z-scores in opposite directions can become highly significant for predictors
that are expected to be correlated in the same direction. We recommend care when
performing the test (i.e. getting many significant genes is not an indication that 
everything worked correctly) and to double-check that marginal TWAS Z-scores are
at least weakly consistent with an association.

---

## FAQ

**How can I validate the TWAS associations?**

We recommend the same procedures for confirming TWAS associations as have been used in GWAS analyses: replication in an external study. This can be done at the level of individual associations, or in aggregate using a gene-based risk score (see manuscript for details). Additionally, if expression weights were computed using this program then the 6th field of the `.imp` file corresponds to the expected accuracy of the prediction on 0-1 scale. Seeing many genes with values well below 0.9 may be an indication of problems with the summary statistics or mismatching LD.

We also proposed a permutation test in the manuscript the evaluates how significant the contribution from expression data is *on top* of the GWAS signal. Loci that pass this permutation test show evidence of heterogeneity that's captured by the expression, and are less likely to be chance co-locoalization.

Lastly, the top eQTL in a gene is expected to explain a substantial fraction of the signal, so individual inspection of eQTLs and conditional analysis can help elucidate the underlying causal variant.

**How should I interpret the effect direction?**

The TWAS effect-size is an estimate of the genetic covariance between gene expression and the GWAS trait, so the direction of effect is informative of this relationship. If the phenotypes have not been treated unusually, then a positive effect means over-expression of the gene leads to an increase in the phenotype (and vice versa). 

**What if the GWAS statistics have highly variable sample sizes?**

In principle having different sample sizes will violate the assumptions made in the IMPG/TWAS approaches and may lead to spurious results. Conservatively, we would recommend to restrict analyses to SNPs that have similar sample sizes, though we have not observed substantial differences in practice. If you do see differences, email us and we can work towards more specific recommendations.

**What else is out there?**

TWAS has analogs to co-localization and Mendelian randomization (with gene expression as one of the traits). [COLOC](https://github.com/chr1swallace/coloc) from the Wallace lab performs Bayesian co-localization analyses. There are many approaches to Mendelian randomization, but [this code](https://github.com/sb452/mr-code) from Stephen Burgess is a great primer; as well as the [SMR/HEIDI](http://cnsgenomics.com/software/smr/) test from the Yang lab. The [MetaXcan](https://github.com/hakyimlab/MetaXcan) and [PrediXcan](https://github.com/hakyimlab/PrediXcan) suite of tools from the Im lab performs gene-based association tests with and without summary data.

You may also be interested in [HESS](http://bogdan.bioinformatics.ucla.edu/software/hess/) to estimate local heritability; [PAINTOR](http://bogdan.bioinformatics.ucla.edu/software/paintor/) to fine-map causal variants; [LD-score Regression](https://github.com/bulik/ldsc) to estimate genome-wide heritability and genetic correlation.

---

## Change Log

| 2016/03/22 | Added omnibus test and weights in `ETC/OMNIBUS`
| 2016/02/24 | Minor update to debug output (monomorphics were not being reported)
| 2016/02/23 | Previous versions had the direction of effect flipped, this did not impact significance but old Z scores must be multipled by -1 to be consistent with individual-level imputation. This has been corrected.
| 2016/02/09 | Updated conversion script to automatically remove monomorphic SNPs
| 2015/12/10 | Updated conversion script to deal with char alleles
| 2015/12/07 | First official release

*Logo by [Ryan Beck](https://thenounproject.com/Ryaaaan/) from The Noun Project*
