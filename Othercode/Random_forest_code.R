# Load necessary libraries
library(randomForest)
library(caret)
library(dplyr)

# Set seed for reproducibility
set.seed(100)

# Split data into training and validation sets
train_indices <- sample(nrow(gene_matrix_for_logi), 0.75 * nrow(gene_matrix_for_logi), replace = FALSE)
TrainSet <- gene_matrix_for_logi[train_indices, ]
ValidSet <- gene_matrix_for_logi[-train_indices, ]

# Train random forest on full feature set
fit_rf <- randomForest(treatment ~ ., data = TrainSet, importance = TRUE)

# View feature importance
RF_importance <- as.data.frame(importance(fit_rf))
RF_importance$gene <- row.names(RF_importance)

# Select top 15 biomarkers based on MeanDecreaseGini
RF_importance_biomarker <- RF_importance %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice(1:15) %>%
  left_join(gene_annotation %>% select(gene, description), by = "gene") %>%
  filter(!is.na(description))

# Plot variable importance
varImpPlot(fit_rf, main = "Top Variable Importance (Gini Index)")

# Prepare reduced training and validation sets
selected_genes <- RF_importance_biomarker$gene

TrainSet_reduced <- TrainSet[, colnames(TrainSet) %in% selected_genes]
TrainSet_reduced$treatment <- TrainSet$treatment

ValidSet_reduced <- ValidSet[, colnames(ValidSet) %in% selected_genes]
ValidSet_reduced$treatment <- ValidSet$treatment

# Train a new random forest model on selected features
model <- randomForest(
  treatment ~ ., 
  data = TrainSet_reduced, 
  ntree = 600, 
  mtry = min(8, length(selected_genes)),  # ensure mtry doesn't exceed number of predictors
  importance = TRUE
)

# Predict on validation set
predValid <- predict(model, ValidSet_reduced, type = "class")

# Plot model error rate
plot(model, main = "Error Rate of Random Forest")

# Evaluate model performance
y_pred <- factor(predValid, levels = c(0, 1))
y_act <- ValidSet_reduced$treatment

# Accuracy
accuracy <- mean(y_pred == y_act)
print(paste("Validation Accuracy:", round(accuracy, 3)))

# Confusion matrix
confusion <- confusionMatrix(y_pred, y_act)
print(confusion)