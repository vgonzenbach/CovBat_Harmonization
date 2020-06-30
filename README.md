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
The R package can be installed via devtools by running the following code

```
# install.packages("devtools")
devtools::install_github("andy1764/CovBat_Harmonization/R")
```

Then, you can load this package via

```
library(CovBat)
```

The R package provides the `covbat` function for harmonization of covariance and the `combat` function which provides some additional options over the original [ComBat](https://github.com/Jfortin1/ComBatHarmonization) package.

For Python, please visit the Python subdirectory.

<div id='id-section2'/>

## 2. Background
Current harmonization methods often focus on addressing scanner differences in the mean and variance of features. However, machine learning methods employed in multivariate pattern analysis (MVPA) are known to leverage additional properties of the data, including covariance. In our recent [preprint](https://www.biorxiv.org/content/10.1101/858415v3), we show that ComBat, a state-of-the-art method designed to harmonize mean and variance, is unable to fully prevent detection of scanner manufacturer through MVPA in the Alzheimer's Disease Neuroimaging Initiative data. We design CovBat to harmonize the covariance of multivariate features and show that it can almost fully prevent detection of scanner properties.

CovBat is meant to be applied after initial preprocessing of the images to obtain a set of features and before statistical analyses. The application of CovBat is not limited to neuroimaging data; however, it has yet to be tested in other types of data.

<div id='id-section3'/>

## 3. Software
The R implementation of CovBat is based on the [ComBat](https://github.com/Jfortin1/ComBatHarmonization) package maintained by Jean-Philippe Fortin. The Python implementation of CovBat is a modification of the ComBat package for Python [here](https://github.com/brentp/combat.py). 