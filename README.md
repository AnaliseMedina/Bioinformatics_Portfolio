# Bioinformatics_Portfolio
This is a collection of the projects I have independently worked on, using the bioinformatic skills I gained from my master's and independent study. 

# Contents 

## Analysis of ALL Dataset in a R markdown

An analysis of expression matrix data that cleans, filters, and normalizes the data to conduct a differentially expressed genes (DEG) analysis. This script utilizes DESeq2 to complete the analysis and other packages for data clustering and visualization. This file should be able to be knit into a html file. 

## Comparison of tidymodels on Breast Cancer Data 

This R script downloads a TCGA breast cancer data set that includes gene expression data and clinical attributes. Samples of interest were chosen before having a few different correlations were determined to see if any genes or attributes correlate with each other. After that, variance analysis and PCA were conducted to further filter the data set before being used to train logistic regression, random forest, SVM, and neural network models whose performances were then compared using AUC and model accuracy. 
