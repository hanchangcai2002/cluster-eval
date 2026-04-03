run_gbmt_analysis <- function(dat, dataset_name, x_vars) {
  results_list <- list()
  dat <- as.data.frame(dat)
  dat$Subject <- as.factor(dat$Subject)
  
  for (k in 2:5) {
    cat("GBMT: Fitting k =", k, "for", dataset_name, "...\n")
    
    start_time <- Sys.time()
    
    fit <- gbmt(
      x.names = x_vars,
      unit = "Subject",
      time = "Time",
      ng = k,
      d = 2,
      data = dat,
      scaling = 0,
      nstart = 3,
      maxit = 5000,
      quiet = FALSE
    )
    
    end_time <- Sys.time()
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
    true_labels <- as.numeric(as.factor(dat$Group[!duplicated(dat$Subject)]))
    pred_labels <- fit$assign
    
    # 计算指标
    ari <- adjustedRandIndex(pred_labels, true_labels)
    bic <- fit$ic["bic"]
    
    # 保存结果
    result <- list(
      dataset = dataset_name,
      method = "gbmt",
      k = k,
      fit = fit,
      metrics = data.frame(
        k = k,
        ARI = ari,
        CHI = NA,
        BIC = bic,
        Runtime_seconds = duration
      )
    )
    results_list[[paste0("k", k)]] <- result
  }
  return(results_list)
}
run_flexmix_analysis <- function(dat, dataset_name, x_vars) {
  results_list <- list()
  model_list <- list()
  for (var in x_vars) {
    model_list[[length(model_list) + 1]] <- FLXMRglm(as.formula(paste(var, "~ Time")))
  }
  dat <- dat[order(dat$Subject, dat$Time), ]
  true_labels <- as.numeric(as.factor(dat$Group[!duplicated(dat$Subject)]))
  
  for (k in 2:5) {
    cat("FlexMix: Fitting k =", k, "for", dataset_name, "...\n")
    
    start_time <- Sys.time()
    
    fit <- flexmix(~1 | Subject, data = dat, k = k, model = model_list)
    
    end_time <- Sys.time()
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
    cluster_subj <- tapply(clusters(fit), dat$Subject, function(x) {
      as.integer(names(which.max(table(x))))
    })
    cluster_subj <- unlist(cluster_subj)
    ari <- adjustedRandIndex(true_labels, cluster_subj)
    bic <- BIC(fit)
    
    result <- list(
      dataset = dataset_name,
      method = "flexmix",
      k = k,
      fit = fit,
      metrics = data.frame(
        k = k,
        ARI = ari,
        CHI = NA,
        BIC = bic,
        Runtime_seconds = duration
      )
    )
    results_list[[paste0("k", k)]] <- result
  }
  
  return(results_list)
}
run_mixAK_analysis <- function(dat, dataset_name, x_vars) {
  results_list <- list()
  dat <- as.data.frame(dat)
  dat <- dat[order(dat$Subject, dat$Time), ]
  y <- dat %>% select(all_of(x_vars)) %>% as.matrix()
  id <- as.numeric(as.factor(dat$Subject))
  
  n_vars <- length(x_vars)
  X_list <- rep(list("empty"), n_vars)
  Z_list <- replicate(n_vars, model.matrix(~ dat$Time - 1), simplify = FALSE)
  randint <- rep(FALSE, n_vars)
  subjects <- unique(dat$Subject)
  n_subjects <- length(subjects)
  n_timepoints <- length(unique(dat$Time))
  true_labels <- as.numeric(as.factor(dat$Group[!duplicated(dat$Subject)]))
  for (k in 2:5) {
    cat("MixAK: Fitting k =", k, "for", dataset_name, "...\n")
    
    start_time <- Sys.time()
    
    fit <- GLMM_MCMC(
      y = y,
      dist = rep("gaussian", n_vars),
      id = id,
      x = X_list,
      z = Z_list,
      random.intercept = randint,
      nMCMC = c(burn = 400, keep = 800, thin = 3, info = 100),
      prior.b = list(Kmax = k),
      init.b = list(K = k),
      PED = TRUE)
    end_time <- Sys.time()
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
    prob_mat <- (fit[[1]]$poster.comp.prob_b + fit[[2]]$poster.comp.prob_b) / 2
    PED <- fit$PED["PED"]
    cluster_pred <- apply(prob_mat, 1, which.max)
    ari <- adjustedRandIndex(cluster_pred, true_labels)
    result <- list(
      dataset = dataset_name,
      method = "mixAK",
      k = k,
      fit = fit,
      cluster_pred = cluster_pred,
      metrics = data.frame(
        k = k,
        ARI = ari,
        CHI = NA,
        PED = PED,
        Runtime_seconds = duration
      )
    )
    results_list[[paste0("k", k)]] <- result
  }
  return(results_list)
}
run_lcmm_analysis <- function(dat, dataset_name, x_vars) {
  results_list <- list()
  dat <- as.data.frame(dat)
  dat$Subject <- as.numeric(dat$Subject)
  formula_str <- paste(paste(x_vars, collapse = " + "), "~ Time")
  cat("LCMM: Initial Model for", dataset_name, "...\n")
  m1 <- multlcmm(
    fixed = as.formula(formula_str),
    random = ~ 1,
    subject = "Subject",
    data = dat,
    link = "linear",
    maxiter = 800,
    verbose = TRUE)
  subjects <- unique(dat$Subject)
  n_subjects <- length(subjects)
  n_timepoints <- length(unique(dat$Time))
  n_vars <- length(x_vars)
  true_labels <- as.numeric(as.factor(dat$Group[!duplicated(dat$Subject)]))
  
  for (k in 2:5) {
    cat("LCMM: Fitting k =", k, "for", dataset_name, "...\n")
    start_time <- Sys.time()
    fit <- multlcmm(
      fixed = as.formula(formula_str),
      mixture = ~ 1,
      random = ~ 1,
      subject = "Subject",
      data = dat,
      link = "linear",
      ng = k,
      B = m1,
      verbose = TRUE,
      maxiter = 800)
    end_time <- Sys.time()
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
    pred_class <- fit$pprob$class
    ari <- tryCatch({
      adjustedRandIndex(pred_class, true_labels)
    }, error = function(e) {
      cat("Warning: ARI calculation failed for k =", k, ":", e$message, "\n")
      return(NA)
    })
    
    bic <- fit$BIC
    ch <- NA
    
    result <- list(
      dataset = dataset_name,
      method = "lcmm",
      k = k,
      fit = fit,
      metrics = data.frame(
        k = k,
        ARI = ari,
        CHI = ch,
        BIC = bic,
        Runtime_seconds = duration
      )
    )
    results_list[[paste0("k", k)]] <- result
  }
  return(results_list)
}

# run_kml3d_analysis <- function(dat, dataset_name, x_vars) {
#   results_list <- list()
#   dat <- as.data.frame(dat)
#   cld3d <- clusterLongData3d(data_array) # 需要相同的时间点, 完全n*t*p的矩阵
#   true_labels <- as.numeric(as.factor(dat$Group[!duplicated(dat$Subject)]))
#   for (k in 2:5) {
#     cat("KML3D: Fitting k =", k, "for", dataset_name, "...\n")
#     start_time <- Sys.time()
#     kml3d(cld3d, nbClusters = k, nbRedrawing = 5, toPlot = "none")
#     end_time <- Sys.time()
#     duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
#     pred_cluster <- getClusters(cld3d, k)
#     ari <- adjustedRandIndex(true_labels, pred_cluster)
#     result <- list(
#       dataset = dataset_name,
#       method = "kml3d",
#       k = k,
#       fit = cld3d,  # 保存完整的clusterLongData3d对象
#       pred_cluster = pred_cluster,  # 单独保存预测结果
#       metrics = data.frame(
#         k = k,
#         ARI = ari,
#         CHI = NA,
#         Runtime_seconds = duration
#       )
#     )
    
#     results_list[[paste0("k", k)]] <- result
#   }
  
#   return(results_list)
# }