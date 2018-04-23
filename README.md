# Multivariate Analysis and Machine Learning on Survival Gene-expression Data

## See "0_Report.pdf" for details.

## Code files
* 1\_Imputation.R: R codes for imputation experiments.   
* 2\_FeatureSelection.ipynb: Python codes for GBDT based feature selection.   
* 3\_Dimensionality\_reduction.R: R codes for the feature extraction experiments.   
* 3.1\_FactorAnalysis.ipynb: Python codes for factor analysis. It is too slow in R, so utilize the Python parallel computation.   
* 3.2\_GaussianMixture.ipynb: Python codes for GMM.   
* 4\_Survival\_Machine\_Learning.ipynb: Python codes for other machine learning parts.   

## Summary

The purpose of this project is to apply various multivariate methods on the survival gene-expression data from [1], and to build better subgrouping and predictor models than the paper. For imputation, Local least squares (LLS) surpassed other 4 models by experiment NRMSE. For dimension reduction, Principal component analysis (PCA) was selected from 16 methods by overall quality scores and visualization. Finally, Gaussian mixture model (GMM) divided patients into 3 subgroups, which has more distinctive Kaplan-Meier curves than the curves of the clustering in [1]. Kernel Survival SVM surpassed other 6 models to achieve the best cross-validation concordance score (0.7398), better than the Cox model (0.6795) in [1].

[1] Rosenwald, Andreas, et al. "The use of molecular profiling to predict survival after chemotherapy for diffuse large-B-cell lymphoma." New England Journal of Medicine 346.25 (2002): 1937-1947.