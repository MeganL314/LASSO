setwd("/Users/lynchmt/Documents/")
## rm(list = ls())
## Load libraries
library(ggplot2)
library(tidyr)
library(rstatix)
library(reshape2)
library(reshape)
library(ggpubr)
library(gridExtra)
library(pROC)
library(nnet)
library(ROCR)
library(survival)
library(survminer)
library(tidyverse)
library(dplyr)
library(plotmo) # for plot_glmnet
### glmnet - Lasso - function cv.glmnet, alpha =1, nfolds=5
library(glmnet)
## https://glmnet.stanford.edu/articles/Coxnet.html#cox-models-for-start-stop-data
## https://stats.stackexchange.com/questions/393937/clarification-for-lasso-based-cox-model-using-glmnet
## lasso <- cv.glmnet(xmat, ysurv, alpha = 1, family = 'cox', nfolds = 5)
## https://www.quantargo.com/help/r/latest/packages/glmnet/makeX.html/cv.glmnet
### https://stats.stackexchange.com/questions/68431/interpretting-lasso-variable-trace-plots
### mlr - Sequential Forward Selection - alpha =0.01, beta=-0.001
library(mlr)
library(MASS)
library(magrittr)
library(cutpointr)
library(survcomp)



## Themes
custom_theme <- function() {
  theme_minimal() %+replace%
    theme(
      
      legend.position = "none",
      plot.title=element_text(hjust=0.5, face="bold", size=16),
      axis.text.y = element_text(size=12),
      axis.text.x = element_text(face="bold", size=12, angle=40),
      axis.title.y = element_text(size=14, angle = 90),
      #margin = margin(t = 0, r = 6, b = 0, l = 0)
      axis.title.x = element_text(size=14, face="bold"),
    )
}

table_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, size=12)
    )
}

custom_bw <- function() {
  theme_bw() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, face="bold", size=18),
      legend.text=element_text(size=12, face="bold"),
      axis.text.y = element_text(size=12),
      axis.text.x = element_text(size=12),
      axis.title.y = element_text(size=14, angle=90),
      axis.title.x = element_text(size=14),
    )
}

probs_theme <- function() {
  theme_minimal() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, size=12),
      legend.position = "none",
      axis.text.y = element_text(size=9, face = "bold"),
      axis.text.x = element_text(size=9, face = "bold"),
      axis.title.y = element_text(size=11, angle = 90, margin = margin(t = 0, r = 4, b = 0, l = 0), face = "bold"),
      axis.title.x = element_text(size=11, face = "bold"),
    )
}

KM_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, face="bold", size=18),
      legend.key.height = unit(.4, 'cm'),
      legend.key.width = unit(.4, 'cm'),
      legend.text=element_text(size=15),
      legend.title=element_text(size=15),
      axis.text.y = element_text(size=16),
      axis.text.x = element_text(size=16),
      axis.title.y = element_text(size=17, angle = 90,  face="bold", margin = margin(t = 0, r = 6, b = 0, l = 0)),
      axis.title.x = element_text(size=17, face="bold"),
      
    )
}


## Helpful functions
create_low_high_list <- function(predicted_data, median) { # create a function with the name my_function
  low_high_list = c()
  for (predicted_prob in predicted_data){
    if (predicted_prob > median){
      low_high_list = c(low_high_list, 'high')
    }
    else if (predicted_prob <= median){
      low_high_list = c(low_high_list, 'low')
    }}
  return(low_high_list)
}

create_labels <- function(median_high, median_low) { # create a function with the name my_function
  if (is.na(median_low)){
    label_low = paste("Low:median = not reached at ", max(dataframe$DFS), " days", sep="")
  }else{label_low = paste("Low:median = ", median_low, " days", sep="")}
  
  if (is.na(median_high)){
    label_high = paste("High:median = not reached at ", max(dataframe$DFS), " days", sep="")
  }else {label_high = paste("High:median = ", median_high, " days", sep="")}
  
  label_list = c(label_high, label_low)
  return(label_list)
}

path <- "QuEST/QuEST_Apply_To_Quicc/"

## Load data
Baseline <- read.csv("./QuEST/OLINK_CBCs_Baseline.csv", header = TRUE, check.names = FALSE)



### Set parameters
dataframe = Baseline
values <- c(1, 0)
index <- c("R", "NR")
dataframe$Responder_num <- values[match(dataframe$Response, index)]

feature_list = c("CD244", "CCL20", "CASP-8")

dataframe[,feature_list] <- lapply(dataframe[,feature_list],as.numeric)


### Set parameters
X_train <- data.matrix(scale(dataframe[,feature_list]))
y_train <- data.matrix(dataframe[,"Responder_num"])


plotname <- 'Olink_limit_3features'

######################################################################################
################################## CHOOSE LAMBDA MIN #################################
######################################################################################

results <- data.frame(seed = integer(), 
                      lambda_min = numeric(), 
                      num_features = integer(), 
                      p_value = numeric(),
                      accuracy_train = numeric(),
                      auc_train = numeric(),
                      'OR>1 Features' = numeric(),
                      'OR<1 Features' = numeric(),
                      selected_analytes = I(list()),
                      check.names=FALSE)  # I(list()) allows list columns



for(seednum in 1:300){
  set.seed(seednum)
  
  # Fit lasso model with cross-validation
  lasso_cv <- cv.glmnet(X_train, y_train, alpha = 1, family = 'binomial', nfolds = 5)
  LambdaMin <- lasso_cv$lambda.min
  
  # Get the coefficients at lambda.min and identify non-zero coefficients
  coef_min <- coef(lasso_cv, s = lasso_cv$lambda.min)
  non_zero_indices <- which(coef_min != 0)  # Indices of non-zero coefficients
  odds_ratios <- exp(coef_min[non_zero_indices])
  
  features_greater_than_1_temp <- ""
  features_less_than_1_temp <- ""
  
  # Flag to control whether to add a comma
  first_greater_than_1 <- TRUE
  first_less_than_1 <- TRUE
  
  for (i in seq_along(non_zero_indices)) {
    feature_name <- rownames(coef_min)[non_zero_indices[i]]
    print(feature_name)
    if (feature_name != "(Intercept)"){
      coef_value <- coef_min[non_zero_indices[i]]
      or_value <- odds_ratios[i]
      print(or_value)
      # Classify features based on OR
      if (or_value > 1) {
        if (first_greater_than_1) {
          features_greater_than_1_temp <- feature_name
          first_greater_than_1 <- FALSE
        } else {
          features_greater_than_1_temp <- paste(features_greater_than_1_temp, feature_name, sep = ", ")
        }}
      
      else if (or_value < 1) {
        print(paste("Odds rato less than 1: ", feature_name, " ", or_value))
        if (first_less_than_1) {
          features_less_than_1_temp <- feature_name
          first_less_than_1 <- FALSE
        } else {
          features_less_than_1_temp <- paste(features_less_than_1_temp, feature_name, sep = ", ")
          
        }
      }}}
  
  num_features <- length(non_zero_indices) - 1  # Exclude intercept

  
  
  # Prepare data for plotting coefficients
  dataframe_coef <- data.frame(coef = coef_min[non_zero_indices], 
                               analyte_full = rownames(coef_min)[non_zero_indices])
  
  # Remove intercept and prepare for ggplot
  dataframe_coef <- dataframe_coef[dataframe_coef$analyte_full != "(Intercept)", ]
  dataframe_coef$s <- ifelse(dataframe_coef$coef < 0, "negative", "positive")
  
  # Sort by absolute value of coefficients (in decreasing order)
  dataframe_coef <- dataframe_coef[order(abs(dataframe_coef$coef), decreasing = TRUE),]
  level_order <- dataframe_coef$analyte_full
  
  # Plot the coefficients using ggplot
  png(paste(path, "/LASSO-out/EffectSizePlots/",plotname, "_", seednum, ".png", sep=""), 
      width=4, height=5, units="in", res=250)
  
  plot <- ggplot(dataframe_coef, aes(x = factor(analyte_full, level = level_order),
                                     y = coef, fill = s)) +
    geom_col(colour = "black") +
    xlab("") +
    ylab("Coefficients") +
    scale_fill_manual(values = c('#B8C0BB', '#D2CBAF'), limits = c("positive", "negative")) +
    ggtitle(paste("lambda = ", round(lasso_cv$lambda.min, 4), sep = "")) +
    custom_theme()
  
  print(plot)
  
  dev.off()
  
  # Get feature names for non-zero coefficients, excluding the intercept
  feature_names <- rownames(coef_min)[non_zero_indices]
  feature_names <- feature_names[feature_names != "(Intercept)"]
  feature_names <- paste(feature_names, collapse = ", ")  # Concatenate feature names with commas

  final_lasso_model <- glmnet(X_train, y_train, alpha = 1, lambda = lasso_cv$lambda.min, family = 'binomial')
  
  # Predict on training and holdout sets using the final fitted model
  preds_train <- predict(final_lasso_model, newx = X_train, type = "response")

  # Convert continuous probabilities to binary predictions (0 or 1) based on 0.5 threshold
  pred_train_binary <- ifelse(preds_train > 0.5, 1, 0)

  # Calculate Accuracy for training 
  accuracy_train <- mean(pred_train_binary == y_train)  # Proportion of correct predictions for training

  # Calculate AUC for training 
  auc_train <- pROC::roc(y_train, preds_train)$auc  # AUC for training set

  # Wilcoxon rank-sum test for p-value
  pval_train <- wilcox.test(preds_train ~ y_train)$p.value


  # Append results
  results <- rbind(results, 
                   data.frame(seed = seednum, 
                              lambda_min = LambdaMin, 
                              num_features = num_features, 
                              p_value = pval_train,
                              accuracy_train = accuracy_train,
                              auc_train = auc_train,
                              'OR>1 Features' = features_greater_than_1_temp,
                              'OR<1 Features' = features_less_than_1_temp,
                              selected_analytes = I(list(feature_names)),
                              check.names=FALSE))
  
}



write.csv(results, paste(path, "/LASSO-out/OR/Olink_QuEST_ForceFeature_LASSO_limit_3features.csv", sep=""), row.names=FALSE)







##### Evaluate the final model
lambda_min = 0.0231 # replace with chosen lambda min

fit <- glmnet(X_train, y_train, family = 'binomial', alpha = 1)
data_coef <- coef(fit, s = lambda_min)
data_coef

####################################################################################
## https://stackoverflow.com/questions/48978179/plotting-lasso-beta-coefficients ##
####################################################################################

# Helper function to plot effect size plots two ways
run_lasso_plots <- function(data_coef, feature_list, input_dataframe, y, path, lambda_min, outfile, survival = FALSE) {
  
  dataframe_coef <- as.data.frame(summary(data_coef))
  
  if (survival == FALSE){
    #remove the intercept column
    row1 <- feature_list[unlist(dataframe_coef$i[2:length(dataframe_coef$i)])-1] 
    row2 <- dataframe_coef$x[2:length(dataframe_coef$x)] #remove the intercept column
    
  } else if (survival == TRUE){
    #remove the intercept column
    row1 <- feature_list[unlist(dataframe_coef$i[1:length(dataframe_coef$i)])] 
    row2 <- dataframe_coef$x[1:length(dataframe_coef$x)] #remove the intercept column
    
  }
  
  
  nonzero_df = data.frame("analyte_full" = row1, "coef" = row2)
  nonzero_df$analyte <- sub("\\.", "-", nonzero_df$analyte_full)
  nonzero_df <- nonzero_df[sort(abs(nonzero_df$coef),decreasing=T,index.return=T)[[2]],]
  nonzero_df$s <- ifelse(nonzero_df$coef < 0, "negative", "positive")
  level_order <- nonzero_df$analyte
  
  important_list = nonzero_df$analyte_full
  input_dataframe[,important_list] <- lapply(input_dataframe[,important_list],as.numeric)
  new_x <- data.matrix(scale(input_dataframe[,important_list]))
  
  if (survival == FALSE){
    fit_elim <- glmnet(new_x, y, family = 'binomial', alpha = 1)
  } else if (survival == TRUE){
    fit_elim <- glmnet(new_x, y, family = 'cox', alpha = 1)
  }
  beta = coef(fit_elim)
  tmp <- as.data.frame(as.matrix(beta))
  tmp$coef <- row.names(tmp)
  tmp <- reshape::melt(tmp, id = "coef")
  tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
  tmp$lambda <- fit_elim$lambda[tmp$variable+1] # extract the lambda values
  tmp$norm <- apply(abs(beta[-1,]), 2, sum)[tmp$variable+1] # compute L1 norm
  tmp$coef = gsub("\\.","-",as.character(tmp$coef))
  
  png(paste(path, outfile,"_plot1.png", sep=""), 
      width=5, height=5, units="in", res=250)
  
  a <- ggplot(tmp[tmp$coef != "(Intercept)",], aes(lambda, value, color = coef, linetype = coef)) + 
    geom_line() + xlab("Lambda") + ylab("Coefficients") +
    guides(color = guide_legend(title = ""), linetype = guide_legend(title = "")) +
    custom_bw() + theme(legend.key.width = unit(3,"lines"))
  print(a)
  dev.off()
  
  png(paste(path, outfile,"_plot2.png", sep=""), 
      width=4, height=5, units="in", res=250)
  
  b <- ggplot(nonzero_df, aes(x=factor(analyte, level = level_order),
                         y=coef, fill=s)) + geom_col(colour="black") + custom_theme() + xlab("") + 
    ylab('Coefficients') + scale_fill_manual(values = c('#B8C0BB', '#D2CBAF'), limits = c("positive", "negative")) +
    ggtitle(paste("lambda = ", round(lambda_min,3), sep=""))
  
  print(b)
  dev.off()
}


run_lasso_plots(data_coef, feature_list, dataframe, y_train, path, lambda_min,
                "LASSO-out/EffectSizePlots/QuEST_Force_3_Feature")

##############################################################
#################### Probability Comparison ###################
##############################################################


lasso_cv <- cv.glmnet(x, y, alpha = 1, family = 'binomial', nfolds = 5)
lasso_cv
preds_1 <- predict(lasso_cv, type="response", family = 'binomial', newx=x, s = lambda_min)
probability_df <- data.frame("Patient" = dataframe[['Sequence Number']], 
                             "Response" = dataframe$Response, 
                             "preds" = c(preds_1))

write.csv(probability_df, 
          paste(path, "/LASSO-out/out-other/QuEST_Force_3_Feature_predprob.csv", sep=""), 
          row.names=FALSE)


png(filename = , width=2.5, height=3.5, units="in", res=250)
my_comparisons <- list(c("R", "NR"))
probability_df$Response <- as.factor(probability_df$Response)
probability_df$Response = factor(probability_df$Response, level=c("NR", "R"))


# Helper function to create plot
plot_probs <- function(df, filename, comparisons, title, ylim_vals) {

  png(filename, width = 2.5, height = 3.5, units = "in", res = 250)
  print(
    ggplot(df, aes(x = Response, y = as.numeric(preds), fill = Response)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = 1.5, height = 0.01, width = 0.2) +
      ylab("Predicted P(Response)") + xlab("") +
      stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                         paired = FALSE, size = 4) +
      ggtitle(title) +
      scale_x_discrete() +
      ylim(ylim_vals) +
      scale_fill_manual(values =c("#A4B279", "#A4B9DA")) +
      probs_theme()
  )
  dev.off()
}

plot_probs(probability_df,
           filename = paste(path, "/LASSO-out/out-other/Olink_QuEST_Force_3_Feature_LASSO_PredProb.png", sep=""),
           comparisons = my_comparisons,
           title = "",
           ylim_vals = c(0, 1.1))



####################################################################
######## Make ROC curve and calculate AUC for LASSO models ########
####################################################################
# library(pROC)
# Assuming `probability_df` contains your predictions and actual responses
roc_curve <- pROC::roc(as.factor(probability_df$Response), probability_df$preds, levels = c("R", "NR"), direction = ">")

# Plot the ROC curve
plot(roc_curve, main="ROC Curve for LASSO Model")

# Calculate AUC
auc_value <- pROC::auc(roc_curve)
print(paste("AUC:", auc_value))

png(paste(path, "/LASSO-out/out-other/QuEST_Force_3_Feature_LASSO_ROCcurve.png", sep=""), width = 5, height = 5, units = "in", res = 300)

plot(roc_curve, 
     main = "ROC Curve for Model Performance", 
     col = "#1c61b6", 
     lwd = 2, 
     xlab = "False Positive Rate", 
     ylab = "True Positive Rate", 
     xlim = c(0, 1), ylim = c(0, 1),
     grid=TRUE)

auc_value <- pROC::auc(roc_curve)
legend("bottomright", legend = paste("AUC =", round(auc_value, 3)), col = "#1c61b6", lwd = 2)

dev.off()





















##########################################################################################################################################################

##########################################################################################################################################################

####################################################################### TEST SET #########################################################################

##########################################################################################################################################################

##########################################################################################################################################################


## Load TEST SET
Quicc <- read.csv("./Quicc/Survival_Models/Quicc_CBC_Olink_Survival_Baseline_Ch28d.csv", header = TRUE, check.names = FALSE)
Arm1 <- Quicc[(Quicc$Arm == "Arm 1"),]

### Set parameters
feature_list = c("CD244_Baseline", "CCL20_Baseline", "CASP.8_Baseline")
Arm1[,feature_list] <- lapply(Arm1[,feature_list],as.numeric)

X_Arm1 <- data.matrix(scale(Arm1[,feature_list]))
y_Arm1 <- Surv(Arm1$survival, as.numeric(as.logical(Arm1$dead)))


plotname <- 'QuiccArm1_limit_features_3Feature_new'

####################################################################################################################################################
########################################################################### CHOOSE LAMBDA MIN ######################################################
#####################################################################################################################################################


results <- data.frame(seed = integer(), 
                      lambda_min = numeric(), 
                      num_features = integer(), 
                      p_value = numeric(),
                      accuracy_train = numeric(),
                      cindex_train = numeric(),
                      'HR>1 Features' = numeric(),
                      'HR<1 Features' = numeric(),
                      selected_analytes = I(list()),
                      check.names=FALSE)  # I(list()) allows list columns



for(seednum in 1:500){
  set.seed(seednum)
  
  # Fit lasso model with cross-validation
  lasso_cv <- cv.glmnet(X_Arm1, y_Arm1, alpha = 1, family = 'cox', nfolds = 5)
  LambdaMin <- lasso_cv$lambda.min
  
  # Get the coefficients at lambda.min and identify non-zero coefficients
  coef_min <- coef(lasso_cv, s = lasso_cv$lambda.min)
  non_zero_indices <- which(coef_min != 0)  # Indices of non-zero coefficients
  odds_ratios <- exp(coef_min[non_zero_indices])  # Hazard ratios (exp of coefficients)
  
  # Temporary variables to store features greater than 1 and less than 1
  features_greater_than_1_temp <- ""
  features_less_than_1_temp <- ""
  
  # Flags to manage adding commas for feature names
  first_greater_than_1 <- TRUE
  first_less_than_1 <- TRUE
  
  # Loop through the non-zero coefficients
  for (i in seq_along(non_zero_indices)) {
    feature_name <- rownames(coef_min)[non_zero_indices[i]]
    if (feature_name != "(Intercept)") {
      coef_value <- coef_min[non_zero_indices[i]]
      hr_value <- odds_ratios[i]
      
      # Print for debugging (optional)
      print(paste("Feature:", feature_name, "HR:", hr_value))
      
      # Classify features based on the HR
      if (hr_value > 1) {
        # If HR > 1, add to the list
        if (first_greater_than_1) {
          features_greater_than_1_temp <- feature_name
          first_greater_than_1 <- FALSE
        } else {
          features_greater_than_1_temp <- paste(features_greater_than_1_temp, feature_name, sep = ", ")
        }
      } else if (hr_value < 1) {
        # If HR < 1, add to the list
        if (first_less_than_1) {
          features_less_than_1_temp <- feature_name
          first_less_than_1 <- FALSE
        } else {
          features_less_than_1_temp <- paste(features_less_than_1_temp, feature_name, sep = ", ")
        }
      }
    }
  }
  
  num_features <- length(non_zero_indices)
  
  
  
  # Prepare data for plotting coefficients
  dataframe_coef <- data.frame(coef = coef_min[non_zero_indices], 
                               analyte_full = rownames(coef_min)[non_zero_indices])
  
  # Remove intercept and prepare for ggplot
  dataframe_coef <- dataframe_coef[dataframe_coef$analyte_full != "(Intercept)", ]
  dataframe_coef$s <- ifelse(dataframe_coef$coef < 0, "negative", "positive")
  
  # Sort by absolute value of coefficients (in decreasing order)
  dataframe_coef <- dataframe_coef[order(abs(dataframe_coef$coef), decreasing = TRUE),]
  level_order <- dataframe_coef$analyte_full
  
  # Plot the coefficients using ggplot
  png(paste(path, "/LASSO-out/EffectSizePlots/",plotname, "_", seednum, ".png", sep=""), 
      width=4, height=5, units="in", res=250)
  
  plot <- ggplot(dataframe_coef, aes(x = factor(analyte_full, level = level_order),
                                     y = coef, fill = s)) +
    geom_col(colour = "black") +
    xlab("") +
    ylab("Coefficients") +
    scale_fill_manual(values = c('#B8C0BB', '#D2CBAF'), limits = c("positive", "negative")) +
    ggtitle(paste("lambda = ", round(lasso_cv$lambda.min, 4), sep = "")) +
    custom_theme()
  
  print(plot)
  
  dev.off()
  
  # Get feature names for non-zero coefficients, excluding the intercept
  feature_names <- rownames(coef_min)[non_zero_indices]
  feature_names <- feature_names[feature_names != "(Intercept)"]
  feature_names <- paste(feature_names, collapse = ", ")  # Concatenate feature names with commas
  
  # **Fit final lasso model using lambda.min**
  final_lasso_model <- glmnet(X_Arm1, y_Arm1, alpha = 1, lambda = lasso_cv$lambda.min, family = 'cox')
  
  # Predict on training and holdout sets using the final fitted model
  preds_train <- predict(final_lasso_model, newx = X_Arm1, type = "response")
  
  # C-index for training set
  cindex_train <- concordance.index(preds_train, surv.time = Arm1$survival, surv.event = as.numeric(as.logical(Arm1$dead)))$c.index


  # **Stratify by median of predicted survival**
  median <- median(preds_train, na.rm = TRUE)
  med_list <- ifelse(preds_train > median, 1, 0)  # Create binary variable for stratification
  med_list <- as.vector(med_list)
  stratify_df <- data.frame("time" = Arm1$survival, "event" = as.numeric(as.logical(Arm1$dead)), "probability" = med_list)
  
  # Fit survival model to stratify based on predicted survival (above and below median)
  fit <- survfit(Surv(time, event) ~ probability, data = stratify_df)
  
  # Calculate p-value for stratified survival curves
  pval_train <- surv_pvalue(fit)$pval
  
  # Append results
  results <- rbind(results, 
                   data.frame(seed = seednum, 
                              lambda_min = LambdaMin, 
                              num_features = num_features, 
                              p_value = pval_train,
                              accuracy_train = accuracy_train,
                              cindex_train = cindex_train,
                              'HR>1 Features' = features_greater_than_1_temp,
                              'HR<1 Features' = features_less_than_1_temp,
                              selected_analytes = I(list(feature_names)),
                              check.names=FALSE))
  
}



write.csv(results, paste(path, "/LASSO-out/HR/Olink_QuiccArm1_ForceFeature_LASSO_3Feature_new.csv", sep=""), row.names=FALSE)







################################################
#### SET SELECTED LAMBDA  ######################
################################################
lambda_min <- 0.0649  # Change to your selected value
plotname <- 'QuiccArm1_3Feature'

################################################
####FIT FINAL MODEL AT SELECTED LAMBDA #########
################################################
fit <- glmnet(X_Arm1, y_Arm1, family = 'cox', alpha = 1, lambda = lambda_min)
coef_min <- coef(fit)
preds_test  <- predict(fit, newx = X_Arm1, type = "response")


run_lasso_plots(coef_min, feature_list, Arm1, y_Arm1, path, lambda_min, outfile, TRUE)

####################################################################################
## https://stackoverflow.com/questions/48978179/plotting-lasso-beta-coefficients ##
####################################################################################

















# KM plots with MEDIAN split

create_labels <- function(data_for_labels, df_followup) {
  label_list <- c()
  group_names <- names(data_for_labels)  # Expected to be something like "group=Lower"
  
  for (i in seq_along(data_for_labels)) {
    group <- sub("^.*=", "", group_names[i])  # Extract just "Lower" or "Upper"
    median_val <- data_for_labels[[i]]
    
    if (is.na(median_val)) {
      label <- paste(group, ": Median = not reached at ", max(df_followup, na.rm = TRUE), " days", sep = "")
    } else {
      label <- paste(group, ": Median = ", median_val, " days", sep = "")
    }
    
    label_list <- c(label_list, label)
  }
  
  return(label_list)
}


plot_km_by_median <- function(preds, survival_time, event_status, filename = NULL) {
  median_val <- median(preds, na.rm = TRUE)
  strat_group <- ifelse(preds > median_val, "Above Median", "Below Median")
  
  # Create binary group by median cutoff
  cutoff <- median(preds, na.rm=TRUE)
  groups <- ifelse(preds >= cutoff, "Upper", "Lower")
  
  # Prepare dataframe
  strat_df <- data.frame(time = survival_time, event = event_status, group = groups)
  colnames(strat_df)[3] <- "group"
  print(strat_df)
  
  # Fit survival curve
  fit <- survfit(Surv(time, event) ~ group, data = strat_df)
  
  # Extract median survival times for each group
  medians <- summary(fit)$table[,"median"]
  
  labels <- create_labels(medians, survival_time)
  
  plot <- ggsurvplot(fit, data = strat_df,
                     pval = TRUE, conf.int = FALSE,
                     risk.table = TRUE,
                     tables.y.text = FALSE,
                     ggtheme = KM_theme(),
                     legend.labs = labels)
  
  plot$plot <- plot$plot + guides(fill = guide_legend(nrow = 2), color = guide_legend(nrow = 2))
  
  
  png(filename, width = 5, height = 5.75, units = "in", res = 300)
  print(plot)
  dev.off()
}

plot_km_by_median(preds_test,
                  Arm1$survival,
                  as.numeric(as.logical(Arm1$dead)),
                  filename = paste0(path, "LASSO-out/out-other/QuiccArm1_3Feature_KM_Plot_Median.png"))











## === C-INDEX: TRAIN + TEST ===

# Train
df_train_cindex <- data.frame(time = y_train_ss$survival,
                              status = as.numeric(as.logical(y_train_ss$dead)),
                              pred = as.numeric(preds_train))

cindex_train_survival <- (1 - concordance(Surv(time, status) ~ pred, data = df_train_cindex)$concordance)


cindex_train_survcomp <- concordance.index(preds_train, surv.time = y_train_ss$survival, 
                                           surv.event = as.numeric(as.logical(y_train_ss$dead)))$c.index

cat("Train C-index (survival::concordance):", round(cindex_train_survival, 2), "\n")
cat("Train C-index (survcomp::concordance.index):", round(cindex_train_survcomp, 2), "\n\n")


# === TEST ===
df_test_cindex <- data.frame(time = y_test_ss$survival, status = as.numeric(as.logical(y_test_ss$dead)),
                             pred = as.numeric(preds_test))

cindex_test_survival <- (1 - concordance(Surv(time, status) ~ pred, data = df_test_cindex)$concordance)

cindex_test_survcomp <- concordance.index(preds_test, surv.time = y_test_ss$survival,
                                          surv.event = as.numeric(as.logical(y_test_ss$dead)))$c.index

cat("Test C-index (survival::concordance):", round(cindex_test_survival, 3), "\n")
cat("Test C-index (survcomp::concordance.index):", round(cindex_test_survcomp, 3), "\n\n")

