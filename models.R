regression <- as.formula ("restrictive ~ age + height + weight + sex +
  fev05 + fev1 + fev3 + fev6 + fvc + fev1_fvc + fef2575 + fefmax +
  expiratory_time"
)

# Logistic regression

model_lr <- train (
  regression,
  data = development,
  method = "glm",
  family = "binomial",
  trControl = trainControl (
    method = "repeatedcv",
    number = 10,
    repeats = 5
  )
)

saveRDS (model_lr, "../models/model_lr.rds")

# Random forest

model_rf <- train (
  regression,
  data = development,
  method = 'ranger',
  importance = 'impurity',
  trControl = trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 5,
    classProbs = T,
    savePredictions = T,
    summaryFunction = twoClassSummary
  ),
  tuneGrid = expand.grid (
    min.node.size = c (20, 30, 40, 50, 60, 70, 80),
    mtry = c (4, 5, 6, 7, 8),
    splitrule = c ("extratrees", "gini", "hellinger")
  )
)

saveRDS (model_rf, "../models/model_rf.rds")

# XGBoost

model_gb <- train (
  regression,
  data = development,
  method = 'xgbTree',
  trControl = trainControl (
    method = "repeatedcv",
    number = 10,
    repeats = 5),
  tuneGrid = expand.grid (
    colsample_bytree = c (0.4, 0.5, 0.6, 0.7, 0.8),
    eta = c (0.001, 0.01, 0.1),
    gamma = c (0.6, 0.7, 0.8, 0.9),
    max_depth = c (2, 4, 6, 8, 10),
    min_child_weight = c (1, 2, 3),
    nrounds = c (400, 600, 800, 1000, 1200),
    subsample = c (0.5, 0.6, 0.7, 0.8, 0.9)
  ),
  verbose = FALSE,
  verbosity = 0
)

saveRDS (model_gb, "../models/model_gb.rds")
