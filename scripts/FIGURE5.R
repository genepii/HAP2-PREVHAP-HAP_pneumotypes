############################################################################################################################
#### CODE USED TO GENERATE FIGURE 5
############################################################################################################################

#A-B-C 

################## Xgboost

library(xgboost)
library(dplyr)
library(readr)
library(tibble)
library(tidyverse)
library(caret)
library(Matrix)
library(pROC)

# Load your data
Bac <- read.delim("bacteriome/virome/combined_fft.txt", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
VAL <- read.delim("validation.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Prepare VAL with outcome placeholder
VAL$C1 <- NA  

set.seed(123)

# Convert train and test C1 to factors and numeric labels
train$C1 <- factor(train$C1)
test$C1 <- factor(test$C1, levels = levels(train$C1)) # maintain level order

train_label <- as.numeric(train$C1) - 1
test_label <- as.numeric(test$C1) - 1

# Select numeric feature columns (exclude C1 and non-numeric)
train_features <- train %>% select(-C1) %>% select(where(is.numeric))
test_features <- test %>% select(-C1) %>% select(where(is.numeric))

# Convert to numeric matrices
train_matrix <- as.matrix(train_features)
test_matrix <- as.matrix(test_features)

# Prepare VAL numeric matrix - keep only features present in train set
VAL_features <- VAL %>% select(colnames(train_features)) %>% select(where(is.numeric))
VAL_matrix <- as.matrix(VAL_features)

# Train XGBoost model
model <- xgboost(
  data = train_matrix,
  label = train_label,
  nrounds = 100,
  objective = "multi:softmax",
  num_class = length(levels(train$C1)),
  verbose = FALSE
)

# Predict on test set
test_pred <- predict(model, test_matrix)

# Calculate and print AUC
roc_obj <- roc(test_label, test_pred)
auc_val <- auc(roc_obj)
cat('AUC:', round(auc_val, 3), '\n')

# Optionally plot ROC curve
plot(roc_obj)

print(table(Predicted = test_pred, Actual = test_label))

cat("\nConfusion Matrix:\n")
conf_mat <- table(Predicted = test_pred, Actual = test_label)
print(conf_mat)

if (length(levels(train$C1)) == 2) {
  # Binary classification metrics
  confusion <- confusionMatrix(
    factor(test_pred, levels = 0:1),
    factor(test_label, levels = 0:1),
    positive = "1"
  )
  
  cat("\nAccuracy:", round(confusion$overall["Accuracy"], 3), "\n")
  cat("Sensitivity:", round(confusion$byClass["Sensitivity"], 3), "\n")
  cat("Specificity:", round(confusion$byClass["Specificity"], 3), "\n")
  
} else {
  # Multi-class accuracy
  accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
  cat("\nMulti-class Accuracy:", round(accuracy, 3), "\n")
}



################## KAPLAN-MEIER AS FIG2A-3I-4I

#D-E --> see FIG3C.R

#F --> see FIG2D-E.R


##############################################################