# CovBat_Harmonization
### Correcting Covariance Batch Effects (CovBat): Harmonization of mean and covariance for multi-site data

--------
**Maintainer**: Andrew Chen, andrewac@pennmedicine.upenn.edu

**License**: Artistic License 2.0

**References**: If you are using CovBat, please cite the ComBat papers:

|       | Citation     | Paper Link
| -------------  | -------------  | -------------  |
| ComBat for multi-site DTI data    | Jean-Philippe Fortin, Drew Parker, Birkan Tunc, Takanori Watanabe, Mark A Elliott, Kosha Ruparel, David R Roalf, Theodore D Satterthwaite, Ruben C Gur, Raquel E Gur, Robert T Schultz, Ragini Verma, Russell T Shinohara. **Harmonization Of Multi-Site Diffusion Tensor Imaging Data**. NeuroImage, 161, 149-170, 2017  |[Link](https://www.sciencedirect.com/science/article/pii/S1053811917306948?via%3Dihub#!)| 
| ComBat for multi-site cortical thickness measurements    | Jean-Philippe Fortin, Nicholas Cullen, Yvette I. Sheline, Warren D. Taylor, Irem Aselcioglu, Philip A. Cook, Phil Adams, Crystal Cooper, Maurizio Fava, Patrick J. McGrath, Melvin McInnis, Mary L. Phillips, Madhukar H. Trivedi, Myrna M. Weissman, Russell T. Shinohara. **Harmonization of cortical thickness measurements across scanners and sites**. NeuroImage, 167, 104-120, 2018  |[Link](https://www.sciencedirect.com/science/article/pii/S105381191730931X)| 
| Original ComBat paper for gene expression array    |  W. Evan Johnson and Cheng Li, **Adjusting batch effects in microarray expression data using empirical Bayes methods**. Biostatistics, 8(1):118-127, 2007.      | [Link](https://academic.oup.com/biostatistics/article/8/1/118/252073/Adjusting-batch-effects-in-microarray-expression) |

If you are using CovBat for harmonization of mean and covariance, please cite the following preprint:

|       | Citation     | Paper Link
| -------------  | -------------  | -------------  |
| Original CovBat paper  | Andrew A. Chen, Joanne C. Beer, Nicholas J. Tustison, Philip A. Cook, Russell T. Shinohara, Haochang Shou, for the Alzheimer’s Disease Neuroimaging Initiative **Removal of Scanner Effects in Covariance Improves Multivariate Pattern Analysis in Neuroimaging Data**. BioRxiv 858415 [Preprint], December 2, 2019. Available from: https://doi.org/10.1101/858415. |[Link](https://www.biorxiv.org/content/10.1101/858415v1)| 

## Table of content
- [1. Installation](#id-section1)
- [2. Background](#id-section2)
- [3. Software](#id-section3)

<div id='id-section1'/>

## 1. Installation
This R package can be installed via devtools by running the following code

```
# install.packages("devtools")
devtools::install_github("andy1764/CovBat_Harmonization")
```

Then, you can load this package via

```
library(CovBat)
```

For Python, please visit the Python subdirectory.

<div id='id-section2'/>

## 2. Background

Placeholder for CovBat introduction. ComBat introduction from https://github.com/Jfortin1/ComBatHarmonization follows:

Imaging data suffer from technical between-scanner variation that hinders comparisons of images across imaging sites, scanners and over time. This includes common imaging modalities, such as MRI, fMRI and DTI, as well as measurements derived from those modalities, for instance ROI volumes, RAVENS maps, cortical thickness measurements, connectome matrices, etc. To maximize statistical power, post-processing data harmonization is a powerful technique to remove unwanted variation when combining data across scanners and sites. 

In two recent papers ([harmonization of DTI data](https://www.sciencedirect.com/science/article/pii/S1053811917306948?via%3Dihub#!) and [harmonization of cortical thickness measurements](https://www.sciencedirect.com/science/article/pii/S105381191730931X)) we have shown that [ComBat](https://academic.oup.com/biostatistics/article/8/1/118/252073/Adjusting-batch-effects-in-microarray-expression), a popular batch-effect correction tool used in genomics, succesffuly removes inter-site technical variability while preserving inter-site biological variability. We showed that ComBat performs well for multi-site imaging studies that only have a few participants per site. We also showed that ComBat was robust to unbalanced studies, that is studies for which the biological covariate of interest is not balanced across sites. 

We recommend to use the ComBat harmonization method after imaging processing, just right before the statistical analysis. The ComBat harmonization requires the imaging data to be represented in a matrix where rows are the imaging features (for instance voxels, ROIs or connectome edges) and columns are the participants. For example, for voxel-level analyses, this usually requires the images to be registered to a common template space. 

The ComBat algorithm needs two mandatory inputs:
- ***The data matrix***. Rows are features and columns are participants. 
- ***The site, study or scanner variable***. The algorithm can only handle one variable. You should provide the smallest unit of the study that you believe introduces unwanted variable. For instance, for a study with 2 sites and 3 scanners (1 site with 1 scanner, 1 site with 2 scanners), the variable for scanner should be used. 

The ComBat algorithm also accepts an optional input:
- ***Biological variables***. You can provide biological covariates, such as disease status, age, gender, to ensure that the harmonization technique does not remove the effects of those variables on the imaging data. The algorithm will take the variability associated with those variables in the estimation of the site/scanner effects. 

<div id='id-section3'/>

## 3. Software

The reference implementation (Standard Version) of ComBat, developed for gene expression analyses, is written in R and is part of the `sva` package available through the Bioconductor project [here](https://bioconductor.org/packages/release/bioc/html/sva.html). This package is an extension of the original ComBat method for harmonization of covariance for multivariate data. We use the same open-source license as the `sva` package, that is the Artistic License 2.0. 