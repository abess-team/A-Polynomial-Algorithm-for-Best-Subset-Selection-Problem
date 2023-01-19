rm(list = ls()); gc(reset = TRUE)
library(snowfall)
source('./source_splicing.R')


## Save the simulation result in file folder: simulation_result
if (!file.exists("simulation_result")) {
  dir.create("simulation_result")
}
setwd("simulation_result")

sfInit(parallel = TRUE, cpus = 6)
sfLibrary(mvtnorm)
sfLibrary(MASS)
sfLibrary(L0Learn)
sfLibrary(glmnet)
sfLibrary(abess)
sfLibrary(BeSS)
sfLibrary(ncvreg)
sfLibrary(EvaluationMeasures)
sfLibrary(microbenchmark)
sfExportAll()

########################################################
###################### low dimension ###################
########################################################
type_seq <- 1
n_seq <- c(40, 60)
sigma_seq <- c(3, 1)
p <- 8
k_seq <- 3
rho_seq <- 0.5
low_dimension_run_time <- FALSE
sfExport("low_dimension_run_time")

M <- 100
for(sigma_value in sigma_seq) {
  for(n in n_seq) {
    res <- sfLapply(1:M, sim_once_low,
                    method = c("abess", "exhaust", "exhaust-sic"),
                    type = 1, n = n, p = p, k = k_seq,
                    rho = rho_seq, sigma = sigma_value, npara = 100)
    res <- simplify2array(res)
    print(paste(sigma_value, n, "done!"))
    res_display <- apply(res, c(1, 2), mean, na.rm = TRUE)
    print(res_display)
    file_name <- paste("low_dimension_original_result", sigma_value, n, sep = "_")
    file_name <- paste0(file_name, ".rda")
    save(res, file = file_name)
  }
}

##########################################################################
################# Runtime Comparison: ABESS v.s. ASR-SIC #################
##########################################################################
type_seq <- 1
n_seq <- c(100)
sigma_seq <- c(1)
p_seq <- 20:40
k_seq <- 3
rho_seq <- 0.5
low_dimension_run_time <- TRUE
sfExport("low_dimension_run_time")

M <- 100
for(p in p_seq) {
  res <- sfLapply(1:M, sim_once_low,
                  method = c("abess", "exhaust", "exhaust-sic"),
                  type = 1, n = n_seq, p = p, k = k_seq,
                  rho = rho_seq, sigma = sigma_seq, npara = 100)
  res <- simplify2array(res)
  res_display <- apply(res, c(1, 2), mean, na.rm = TRUE)
  print(paste(p, "done!"))
  print(res_display)
  file_name <- paste("runtime_original_result", p, sep = "_")
  file_name <- paste0(file_name, ".rda")
  save(res, file = file_name)
}

########################################################
##################### high dimension ###################
########################################################
low_dimension_run_time <- FALSE
sfExport("low_dimension_run_time")
SNR_seq <- 5
type_seq <- c(1, 2)
n_seq <- c(500)
p_seq <- c(500, 1500, 2500)
k_seq <- 10
rho_seq <- 0.8

M <- 200
for(type in type_seq) {
  for(p in p_seq) {
    para_num <- floor(n_seq / (log(p) * log(log(n_seq))))
    res <- sfLapply(1:M, sim_once_high, SNR = SNR_seq, type = type,
                    n = n_seq, p = p, k = k_seq,
                    rho = rho_seq, sigma = 1, npara = para_num)
    res <- simplify2array(res)
    print(paste(type, p, "done!"))
    file_name <- paste("high_dimension_original_result", type, p, sep = "_")
    file_name <- paste0(file_name, ".rda")
    save(res, file = file_name)
  }
}
sfStop()
