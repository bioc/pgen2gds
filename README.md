pgen2gds: Format Conversion between PLINK2 PGEN and GDS
===

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)


## Description

This package provides functions for format conversion from [PLINK2 pgen](https://www.cog-genomics.org/plink/2.0) files to [SeqArray GDS](https://www.bioconductor.org/packages/SeqArray) files, and from SeqArray GDS files to PLINK2 pgen files.


## Version

v0.99.4


## Package Maintainer

Dr. Xiuwen Zheng


## Installation

Requires R (≥ v4.0.0), [gdsfmt](http://www.bioconductor.org/packages/gdsfmt), [SeqArray](http://www.bioconductor.org/packages/SeqArray)

* Installation from Biocondutor:
```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("pgen2gds")
```

* Installation from Github:
```R
library(remotes)
install_github("zhengxwen/pgen2gds")
```
The `install_github()` approach requires that you build from source, i.e. `make` and compilers must be installed on your system -- see the [R FAQ](http://cran.r-project.org/faqs.html) for your operating system; you may also need to install dependencies manually.


## Citations for GDS

Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS (2012). A High-performance Computing Toolset for Relatedness and Principal Component Analysis of SNP Data. *Bioinformatics*. [DOI: 10.1093/bioinformatics/bts606](http://dx.doi.org/10.1093/bioinformatics/bts606).

Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance data format for WGS variant calls. *Bioinformatics*. [DOI: 10.1093/bioinformatics/btx145](http://dx.doi.org/10.1093/bioinformatics/btx145).


## Examples

```R
library(pgen2gds)

# Example file
pgen_fn <- system.file("extdata", "plink2_gen.pgen", package="pgen2gds")

# Format conversion
seqPGEN2GDS(pgen_fn, out.gdsfn="test.gds")
## PLINK2 PGEN to SeqArray GDS:
##     pgen file (14.3K):
##         plink2_gen.pgen
##     pvar file (12.7K):
##         plink2_gen.pvar
##     psam file (34.3K):
##         plink2_gen.psam
##         reading ...
##     # of samples: 2504
##     # of variants: 482
##     Output:
##         test.gds
## ...

# Format conversion from GDS to PGEN (pgen, pvar and psam files)
seqGDS2PGEN("test.gds", "test_out")
## SeqArray GDS to PLINK2 PGEN:
##     open 'test.gds'
##     # of samples: 2,504
##     # of variants: 482
##     [Output]
##     PSAM: test_out.psam
##     PVAR: test_out.pvar
##     PGEN: test_out.pgen
## ...
```

The pgen writer stores each variant record in the smallest of the dense 2-bit, 1-bit and difference-list encodings (with multiallelic and phase tracks), following the thresholds used by PLINK2, but it does not use the LD-based encodings of PLINK2. The output can therefore be slightly larger than a pgen file written by PLINK2 itself; run `plink2 --pfile test_out --make-pgen --out new` if the smallest file size is needed.


## Also See

[seqVCF2GDS()](https://rdrr.io/bioc/SeqArray/man/seqVCF2GDS.html) in the [SeqArray](https://bioconductor.org/packages/SeqArray) package, conversion from VCF files to GDS files.

[seqBED2GDS()](https://rdrr.io/bioc/SeqArray/man/seqBED2GDS.html) in the [SeqArray](https://bioconductor.org/packages/SeqArray) package, conversion from PLINK BED files to GDS files.

[seqBGEN2GDS()](https://github.com/CoreArray/gds2bgen) conversion from BGEN files to GDS files.
