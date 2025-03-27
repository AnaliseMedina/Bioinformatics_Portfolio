# Loading in required libraries 

library(TCGAbiolinks); library(EDASeq)
library(SummarizedExperiment)
library(dplyr)
library(tidyverse)
library(reshape2)
library(corrplot)
library(factoextra)
library(tidyr)
library(tidymodels)
library(rsample)
library(vip)
library(randomForest)
library(future)
library(ranger)
library(AppliedPredictiveModeling)
library(brulee)
library(torch)
library(tibble)
library(yardstick)

# Downloading TCGA dataset

## Obtaining Project Summary
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")

## Downloading and preparing Target Case

### Querying available samples
TargetSamples <- GDCquery(project = "TCGA-BRCA", 
                          data.category = "Transcriptome Profiling", 
                          data.type = "Gene Expression Quantification",
                          workflow.type = "STAR - Counts")

CaseInfo <- getResults(TargetSamples)

### Subsetting the data so there is an equal amount of control and cancerous samples 
dataPrimary_Target <- TCGAquery_SampleTypes(barcode = CaseInfo$cases, typesample = "TP") # primary tumor
dataNormal_Target <- TCGAquery_SampleTypes(barcode = CaseInfo$cases, typesample = "NT") # normal tissue
dataPrimary_Target <- dataPrimary_Target[1:113]
dataNormal_Target <- dataNormal_Target[1:113]


## Downloading the samples of interest
TargetSamples <- GDCquery(project = "TCGA-BRCA",
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification",
                          workflow.type = "STAR - Counts",
                          barcode = c(dataPrimary_Target, dataNormal_Target))

GDCdownload(TargetSamples)


#Preparing the data 

data <- GDCprepare(TargetSamples)

## Filtering for protein-coding genes
table(rowData(data)$gene_type)

SECoding <- data[rowData(data)$gene_type == "protein_coding", ]

## Removing low-quality samples 
dataPrep_Coding <- TCGAanalyze_Preprocessing(object = SECoding, cor.cut = 0.6,  datatype = "fpkm_unstrand")

boxplot(dataPrep_Coding, outline = FALSE)

## Normalizing and further filtering the data
dataNorm_Coding <- TCGAanalyze_Normalization(tabDF = dataPrep_Coding, geneInfo = geneInfoHT, method = "geneLength")
dataFilt_Coding <- TCGAanalyze_Filtering(tabDF = dataPrep_Coding, method = "quantile", qnt.cut =  0.25)   
boxplot(dataNorm_Coding, outline = FALSE)

## Differential Expression Analysis 
DEGsCoding <- TCGAanalyze_DEA(mat1 = dataFilt_Coding[,dataNormal_Target],
                              mat2 = dataFilt_Coding[,dataPrimary_Target],
                              pipeline="limma",
                              Cond1type = "Normal",
                              Cond2type = "Tumor",
                              fdr.cut = 0.01 ,
                              logFC.cut = 1,
                              method = "glmLRT", ClinicalDF = data.frame())

# Correlating genes to continuous clinical attributes 

## Age and gene expression 

clin_data <- data.frame(data@colData) 
sample_ids <- colnames(dataFilt_Coding)

clin_data <- clin_data[rownames(clin_data) %in% sample_ids, ] #Filtering through the clinical data

age_data <- clin_data[,"age_at_diagnosis", drop = FALSE] #Obtaining only the age data and turning it into years
age_data <- round(age_data/365)
age_data <- na.omit(age_data)

t_ex <- as.data.frame(t(dataFilt_Coding)) #Transposing so the dimensions match that of the age_data 
t_ex_1 <- t_ex[rownames(t_ex) %in% rownames(age_data),]

#Correlating Genes to Continuous Attributes 

## Genes vs. Age at Diagnosis

age_corr <- numeric(ncol(t_ex_1)) #Place to put correlation data

for (i in 1:ncol(t_ex_1)){
  age_corr[i] <- cor(t_ex_1[,i], age_data, use = "complete.obs")
}

age_corr_df <- data.frame(gene = rownames(t(t_ex_1)),
                          correlation = age_corr)

age_corr_df <- age_corr_df[order(age_corr_df$correlation, decreasing = FALSE),] #Negative Correlations 
head(age_corr_df)

age_corr_df <- age_corr_df[order(age_corr_df$correlation, decreasing = TRUE),] #Positive Correlations
head(age_corr_df)

hist(age_corr_df$correlation, main = "Correlation between Gene Expression and Age at Diagnosis", xlab = "Correlation Values")

## Initial weight and gene expression

dtd_data <- clin_data[,"initial_weight", drop = FALSE] 
dtd_data <- na.omit(dtd_data)

t_ex_2 <- t_ex[rownames(t_ex) %in% rownames(dtd_data),]

dtd_corr <- numeric(ncol(t_ex_2)) #Place to put correlation data

for (i in 1:ncol(t_ex_2)){
  dtd_corr[i] <- cor(t_ex_2[,i], dtd_data, use = "complete.obs")
}


dtd_corr_df <- data.frame(gene = rownames(t(t_ex_2)),
                          correlation = dtd_corr)

dtd_corr_df <- dtd_corr_df[order(dtd_corr_df$correlation, decreasing = FALSE),] #Negative Correlations 
head(dtd_corr_df)

dtd_corr_df <- dtd_corr_df[order(dtd_corr_df$correlation, decreasing = TRUE),] #Positive Correlations
head(dtd_corr_df)

hist(dtd_corr_df$correlation, main = "Correlation between Gene Expression and Initial Weight", xlab = "Correlation Values")

## Pathological stage and gene expression

clin_data <- clin_data %>%
  mutate(paper_pathologic_stage = case_when(
    paper_pathologic_stage == "Stage_I" ~ "1",
    paper_pathologic_stage == "Stage_II" ~ "2", 
    paper_pathologic_stage == "Stage_III" ~ "3",
    paper_pathologic_stage == "Stage_IV" ~ "4",
    TRUE ~ paper_pathologic_stage
  ))

stage_data <- clin_data[,"paper_pathologic_stage", drop = FALSE] 
stage_data <- na.omit(stage_data)
stage_data$paper_pathologic_stage <- as.numeric(stage_data$paper_pathologic_stage)


t_ex_3 <- t_ex[rownames(t_ex) %in% rownames(stage_data),]


stage_corr <- numeric(ncol(t_ex_3)) #Place to put correlation data

for (i in 1:ncol(t_ex_3)){
  stage_corr[i] <- cor(t_ex_3[,i], stage_data, use = "complete.obs")
}


stage_corr_df <- data.frame(gene = rownames(t(t_ex_3)),
                            correlation = stage_corr)

stage_corr_df <- stage_corr_df[order(stage_corr_df$correlation, decreasing = FALSE),] #Negative Correlations 
head(stage_corr_df)

stage_corr_df <- stage_corr_df[order(stage_corr_df$correlation, decreasing = TRUE),] #Positive Correlations
head(stage_corr_df)

hist(stage_corr_df$correlation, main = "Correlation between Gene Expression and Pathological Stage", xlab = "Correlation Values")

# Correlating gene expression between genes 

cor_matrix <- cor(t_ex, method = "pearson")
cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- NA # Removing duplicates 

cor_long <- melt(cor_matrix, na.rm=TRUE)
top_cor <- cor_long[order(-cor_long$value),][1:1000,]
head(top_cor) # Most Correlated Genes 

top_cor_genes <- unique(c(top_cor$Var1, top_cor$Var2)) 

hist(cor_long$value) #Observing distribution of correlation coefficients
summary(cor_long$value)

count <- sum(cor_long$value >= 0.9); count #Seeing how many genes have extremely high correlation coefficients

gene_count <- table(top_cor$Var1)
head(sort(gene_count, decreasing = TRUE))

gene_count <- table(top_cor$Var2)
head(sort(gene_count, decreasing = TRUE))

# Correlating malignancy, pathological stage, age at diagnosis, and initial weight

meta_corrs <- data.frame(
  malignancy = clin_data$synchronous_malignancy,
  stage = clin_data$paper_pathologic_stage,
  age = round(clin_data$age_at_diagnosis/365),
  weight = clin_data$initial_weight)

rownames(meta_corrs) <- rownames(clin_data)

meta_corrs <- na.omit(meta_corrs)

meta_corrs <- meta_corrs[rownames(meta_corrs ) != "TCGA-A8-A099-01A-11R-A00Z-07", ] #Removing this specific row because the NA would not be removed by the previous line. 

meta_corrs$malignancy <- factor(meta_corrs$malignancy, levels = c("No", "Yes"))

meta_corrs$stage <- as.numeric(meta_corrs$stage)

meta_corrs$malignancy <- as.numeric(meta_corrs$malignancy) - 1

meta_correlations <- cor(meta_corrs, method = "spearman") #Used this instead of Pearson in order to be less sensitive to outliers 

corrplot(meta_correlations, method = "color") #Constructing Correlogram

# Variance Analysis 

boxplot(dataFilt_Coding, main = "Expression Variance Across Samples", las = 2)

gene_var <- apply(dataFilt_Coding,1,var)
gene_means <- apply(dataFilt_Coding, 1, mean)

plot(gene_means, gene_var, log="xy", pch=16, col="blue", 
     xlab = "Mean Expression", ylab = "Variance",
     main = "Mean-Variance Plot")
abline(a=0,b=1, col = "purple", lwd = 2)

top_genes <- names(sort(gene_var, decreasing = TRUE)[1:2000])
head(top_genes) # Genes with the greatest variance

bottom_genes <- names(sort(gene_var, decreasing = FALSE))
head(bottom_genes) # Genes with the least variance 

ex <- dataFilt_Coding[top_genes,] #To use for all the upcoming PCA. 

t_ex2 <- t(ex)

# PCA analysis

new_t_ex <- t(ex)

scaled_t_ex <- scale(new_t_ex)

pca_res <- prcomp(scaled_t_ex, scale. = TRUE)

fviz_pca_ind(pca_res,
             col.ind = clin_data$tissue_type,  
             palette = "jco",            
             title = "PCA - Tumor vs. Normal",
             label = "none") +          
  scale_color_manual(values = c("red", "blue"))

fviz_pca_ind(pca_res,
             col.ind = clin_data$paper_pathologic_stage,  
             palette = "jco",            
             title = "PCA - Pathological Stage",
             label = "none") +          
  scale_color_manual(values = c("red", "blue", "green", "purple", "grey"))

fviz_pca_ind(pca_res,
             col.ind = clin_data$synchronous_malignancy,  
             palette = "jco",            
             title = "PCA - Malignancy",
             label = "none") +          
  scale_color_manual(values = c("red", "blue")) 

fviz_pca_ind(pca_res,
             col.ind = clin_data$paper_BRCA_Subtype_PAM50,  
             palette = "jco",            
             title = "PCA - Subtype",
             label = "none") +          
  scale_color_manual(values = c("red", "blue","purple", "darkslateblue", "green", "grey")) 

pca_loading <- pca_res$rotation

top_genes <- names(sort(abs(pca_loading[,1]), decreasing = TRUE))[1:1000]

ex <- dataFilt_Coding[top_genes,]
new_t_ex <- t(ex)

gc()

# Model Development and comparison

## Data splitting 

groups <- clin_data$tissue_type

data <- cbind(groups, new_t_ex) %>% as.data.frame() %>% mutate(groups = factor(groups, levels = c("Normal", "Tumor")))

data <- data %>%
  mutate(across(-groups, as.numeric))  # Convert all columns except 'groups' to numeric

data_split <- initial_split(data, strata = "groups")
data_train <- training(data_split)
data_test <- testing(data_split)

folds <- vfold_cv(data_train, strata = groups, repeats = 10) #Kept the repeats low in interest of computational time, but could have increased folds and repeats for statistical vigor. 

## Creating Recipe

recipe <- recipe(data_train) %>%
  update_role(colnames(data), new_role="predictor")  %>%
  update_role(groups, new_role="outcome") %>%
  step_zv(all_numeric(), -all_outcomes()) %>%
  step_normalize(all_numeric(), -all_outcomes())

## Metrics Dataframe

accuracy_df <- tibble(
  model = character(),
  accuracy = numeric(),
  roc_auc = numeric()
)

## Logistic Regression 

logreg_spec <- logistic_reg(penalty = tune(), mixture = tune()) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

logreg_workflow <- workflow() %>%
  add_recipe(recipe) %>%
  add_model(logreg_spec)

log_reg_grid <- grid_regular(
  penalty(range = c(-5, 0), trans = log10_trans()), # Log scale for better tuning
  mixture(range = c(0, 1)), # Explore Ridge to LASSO
  levels = 10
)

#run the tuning grid
logreg_grid <- tune_grid(
  logreg_workflow,
  resamples = folds,
  grid=log_reg_grid
)

#get the tuning metrics
logreg_metrics<-logreg_grid %>%
  collect_metrics()

highest_roc_auc <- logreg_grid %>%
  select_best(metric="accuracy")  

final_log_reg <- finalize_workflow(logreg_workflow, highest_roc_auc)

logreg_fit <- fit(final_log_reg, data = data_train)

#make class predictions
class_preds <- predict(logreg_fit, new_data = data_test,
                       type = 'class')
#make probability predictions
prob_preds <- predict(logreg_fit, new_data = data_test,
                      type = 'prob')

#combine test y and results into a dataframe
logreg_results<- data_test %>%
  select(groups) %>%
  bind_cols(class_preds, prob_preds)

#calculate the AUC
auc<-roc_auc(logreg_results,
             truth = groups,
             ".pred_Normal",
             event_level="first")$.estimate

#confustion matrix
conf_mat(logreg_results, truth = groups, estimate = ".pred_class")

#get classification metrics
classification_metrics <- metric_set(accuracy, f_meas, spec, sens, npv, ppv)
classification_metrics(logreg_results, truth = groups, estimate = ".pred_class")

acc_value <- logreg_results %>%
  accuracy(truth = groups, estimate = .pred_class) %>%
  pull(.estimate)

accuracy_df <- add_row(accuracy_df, 
                       model = "Logistic Regression",
                       accuracy = acc_value,
                       roc_auc = auc)

#generate an ROC curve
g_roc<- logreg_results %>%
  roc_curve(truth=groups, paste0(".pred_Normal"), event_level="first") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path(color="red") +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw() +
  annotate(geom="text", x=0.75, y=0.1, label=paste0("AUC ", round(auc, 3)), color="red") +
  ggtitle(paste0("logreg ROC")); g_roc

g_importance_vip<-logreg_fit%>%
  extract_fit_parsnip() %>%
  vip::vip(geom = "col", num_features = 20, mapping = aes(fill = Sign)) +
  theme_bw(base_size = 14) +
  labs(title = paste0("logreg Importance")) 

## Random Forest

tune_spec <- rand_forest(trees = tune(), mtry = tune(), min_n = tune()) %>%
  set_engine("ranger", importance = "impurity") %>% 
  set_mode("classification")

#mtry is dependent on number of columns
rf_param<-extract_parameter_set_dials(tune_spec) %>%
  finalize(data_train) 

#set a tree grid
tree_grid<-grid_regular(mtry(range = c(2, 10)),  
                        trees(range = c(100, 500)), 
                        min_n(range = c(2, 10)),   
                        levels = 5                
)



#update the workflow spec with the new model
tune_workflow <- workflow() %>%
  add_recipe(recipe) %>%
  add_model(tune_spec) 

set.seed(1234)
#system start time
a<-Sys.time()

plan(multisession, workers = 4)  


#run the tuning grid
rf_grid <- tune_grid(
  tune_workflow,
  resamples = folds,
  grid=tree_grid
)
#system end time
b<-Sys.time()
b-a

#get the tuning metrics
rf_metrics<-rf_grid %>%
  collect_metrics()

highest_roc_auc <- rf_grid %>%
  select_best(metric="accuracy")  

#finalize the model
final_rf <- finalize_model(
  tune_spec,
  highest_roc_auc
)

#finalize the workflow
final_wf <- workflow() %>%
  add_recipe(recipe) %>%
  add_model(final_rf)

#fit the final tuned model
rf_tune_fit <- fit(final_wf, data = data_train)

#predict class, probability
class_preds_tune <- predict(rf_tune_fit, new_data = data_test,
                            type = 'class')
prob_preds_tune <- predict(rf_tune_fit, new_data = data_test,
                           type = 'prob')

#collate results
rf_tune_results<- data_test %>%
  select(groups) %>%
  bind_cols(class_preds_tune, prob_preds_tune)

#calculate AUC
auc_tune<-roc_auc(rf_tune_results,
                  truth = groups,
                  ".pred_Normal", 
                  event_level="first")$.estimate

acc_value <- rf_tune_results %>%
  accuracy(truth = groups, estimate = .pred_class) %>%
  pull(.estimate)

accuracy_df <- add_row(accuracy_df, 
                       model = "Random Forest",
                       accuracy = acc_value,
                       roc_auc = auc_tune)

#confusion matrix and metrics
conf_mat(rf_tune_results, truth = groups, estimate = ".pred_class")
classification_metrics <- metric_set(accuracy, f_meas, spec, sens, npv, ppv)
classification_metrics(rf_tune_results, truth = groups, estimate = ".pred_class")


#tune ROC curve
g_roc_tune<- rf_tune_results %>%
  roc_curve(truth=groups, paste0(".pred_Normal"), event_level="first") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path(color="red") +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw() +
  annotate(geom="text", x=0.75, y=0.1, label=paste0("AUC ", round(auc_tune, 3)), color="red") +
  ggtitle(paste0("Tuned rf ROC")); g_roc_tune

## SVM 

svm_grid <- grid_regular(
  cost(range = c(-5, 5)),  # Log scale for C
  rbf_sigma(range = c(-5, 2)),  # For RBF kernel
  levels = 10
)


#in the model spec, change cost to tune()
tune_spec <- svm_rbf(cost = tune(), rbf_sigma=tune()) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

#update the workflow spec with the new model
tune_workflow <- workflow() %>%
  add_recipe(recipe) %>%
  add_model(tune_spec) 

set.seed(1234)
#system start time
a<-Sys.time()
#run the tuning grid
svm_grid <- tune_grid(
  tune_workflow,
  resamples = folds,  
  grid=svm_grid
)
#system end time
b<-Sys.time()
b-a

#get the tuning metrics
svm_metrics<-svm_grid %>%
  collect_metrics()

#pick the best cost on highest roc_auc or accuracy
highest_roc_auc <- svm_grid %>%
  select_best(metric="accuracy")  

#finalize the model
final_svm <- finalize_model(
  tune_spec,
  highest_roc_auc
)

#finalize the workflow
final_wf <- workflow() %>%
  add_recipe(recipe) %>%
  add_model(final_svm)

#fit the final tuned model
svm_tune_fit <- fit(final_wf, data = data_train)

#predict class, probability
class_preds_tune <- predict(svm_tune_fit, new_data = data_test,
                            type = 'class')
prob_preds_tune <- predict(svm_tune_fit, new_data = data_test,
                           type = 'prob')

#collate results
svm_tune_results<- data_test %>%
  select(groups) %>%
  bind_cols(class_preds_tune, prob_preds_tune)

#calculate AUC
auc_tune<-roc_auc(svm_tune_results,
                  truth = groups,
                  ".pred_Normal", 
                  event_level="first")$.estimate

acc_value <- svm_tune_results %>%
  accuracy(truth = groups, estimate = .pred_class) %>%
  pull(.estimate)

accuracy_df <- add_row(accuracy_df, 
                       model = "SVM",
                       accuracy = acc_value,
                       roc_auc = auc_tune)

#confusion matrix and metrics
conf_mat(svm_tune_results, truth = groups, estimate = ".pred_class")
classification_metrics <- metric_set(accuracy, f_meas, spec, sens, npv, ppv)
classification_metrics(svm_tune_results, truth = groups, estimate = ".pred_class")


#tune ROC curve
g_roc_tune<- svm_tune_results %>%
  roc_curve(truth=groups, paste0(".pred_Normal"), event_level="first") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path(color="red") +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw() +
  annotate(geom="text", x=0.75, y=0.1, label=paste0("AUC ", round(auc_tune, 3)), color="red") +
  ggtitle(paste0("Tuned SVM ROC")); g_roc_tune

## Neural Network 

nnet_spec<- mlp(epochs = 1000, hidden_units = c(8,16), 
                penalty = 0.01, learn_rate = 0.1) %>%
  set_engine("brulee", importance = "permutation") %>% 
  set_mode("classification")

#create the workflow from the recipe and the model
nnet_workflow <- workflow() %>%
  add_recipe(recipe) %>%
  add_model(nnet_spec)

set.seed(31416)
#system start time
a<-Sys.time()
#fit the training data to the workflow
nnet_fit <- fit(nnet_workflow, data = data_train)
#system end time
b<-Sys.time()
#evaluate time
b-a

#make class predictions
class_preds <- predict(nnet_fit, new_data = data_test,
                       type = 'class')
#make probability predictions
prob_preds <- predict(nnet_fit, new_data = data_test,
                      type = 'prob')

#combine test y and results into a dataframe
nnet_results<- data_test %>%
  select(groups) %>%
  bind_cols(class_preds, prob_preds)

#calculate the AUC
auc<-roc_auc(nnet_results,
             truth = groups,
             ".pred_Normal", #ground truth
             event_level="first")$.estimate

acc_value <- nnet_results %>%
  accuracy(truth = groups, estimate = .pred_class) %>%
  pull(.estimate)

accuracy_df <- add_row(accuracy_df, 
                       model = "Neural Network",
                       accuracy = acc_value,
                       roc_auc = auc)

#confustion matrix
conf_mat(nnet_results, truth = groups, estimate = ".pred_class")

#get classification metrics
classification_metrics <- metric_set(accuracy, f_meas, spec, sens, npv, ppv)
classification_metrics(nnet_results, truth = groups, estimate = ".pred_class")

#generate an ROC curve
g_roc<- nnet_results %>%
  roc_curve(truth=groups, paste0(".pred_Normal"), event_level="first") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path(color="red") +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw() +
  annotate(geom="text", x=0.75, y=0.1, label=paste0("AUC ", round(auc, 3)), color="red") +
  ggtitle(paste0("nnet ROC")); g_roc

## Comparing accuracy and AUC between these different models 

accuracy_df
