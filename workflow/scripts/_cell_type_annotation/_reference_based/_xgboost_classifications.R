#' Run xgboost regression for a set of genes
#'
#' Train model for cell-type classification based on panel genes
#'
#' @param X reference (i.e., gene expression matrix or low-dim representation of genes to model from) in a format **cells x features**
#' @param Y label (i.e, cell types) in a format **targets x cells**
#' @param targets a vector of target genes
#' @param xgb_objective \code{objective} of \link[xgboost]{xgb.train} function
#' @param do.cross.val whether to run cross validation
#' @param train.index cell indices to train model for cross validation (random half if not provided)
#' @param num_class number of classes (cell types)
#' @param ... other parameters from \link[xgboost]{xgb.train}
#'
#' @return list of \link{run_regression_per_gene} results with names being \code{targets}
#'

run_regressions_classification_fixed_features <- function(
    X,
    Y,
    do.cross.val = TRUE,
    train.index = NULL,
    xgb_objective = "multi:softprob",
    nrounds = 1000,
    max_depth = 2,
    eta = 0.03,
    nthread = 1,
    verbose = FALSE, 
    num_class = NULL,
    ...
) {
  res <- list()
  if(is.null(num_class))
    num_class <- length(unique(Y))
  
  if (do.cross.val & is.null(train.index)) {
    warning("`train.index` not provided, setting random half...")
    set.seed(123)
    cell.ids      <- rownames(X)
    N.train       <- round(length(cell.ids) * 0.5 + 1)
    train.index   <- sample(length(cell.ids), N.train, replace = F)
  }
    
  if (do.cross.val) {
    ## checking to add
    data.train <- xgb.DMatrix(data = X[train.index,],
                              label = Y[train.index])
    
    data.test  <- xgb.DMatrix(data = X[-train.index,],
                              label = Y[-train.index])
    watchlist <- list(train = data.train, test = data.test)
    
    set.seed(1) # needed?
    fit.train <- xgb.train(
      data = data.train,
      max_depth = max_depth,
      watchlist = watchlist,
      eta = eta,
      nrounds = nrounds,
      objective = xgb_objective,
      verbose = verbose,
      num_class = num_class,
      ...
    )
   # nrounds <-
   #   which.min(fit.train$evaluation_log[[3]])
    
    res[["fit.train"]] <- fit.train
  }
  
  fit <- xgboost(
    data = X,
    label = Y,
    objective = xgb_objective,
    max_depth = max_depth,
    eta = eta,
    nrounds = nrounds,
    nthread = nthread,
    verbose = verbose,
    num_class = num_class,
    ...
  )
  
  res[["fit"]] <- fit
  
  return(res)
}