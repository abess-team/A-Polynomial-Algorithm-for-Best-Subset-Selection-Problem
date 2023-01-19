rm(list = ls()); gc(reset = TRUE)
source("boldlatex.R")
library(reshape2)
library(dplyr)
library(abind)
library(xtable)

setwd("./simulation_result/")

######################################################################
############################## Table 1 ###############################
######################################################################
pre_process_low <- function(res) {
  selected_metric <- c(5, 6, 8, 4)
  selected_method <- 1:3
  res_tmp <- res[, selected_metric, ]
  res_tmp <- res_tmp[selected_method, , ]
  res_tmp[, 2, ] <- 1 - res_tmp[, 2, ]
  res_tmp[, 3, ] <- (res_tmp[, 3, ] - 3)
  res_tmp <- res_tmp[, c(1, 2, 4, 3), ]
  res_tmp
}

load("low_dimension_original_result_3_40.rda")
processed_res_case1 <- pre_process_low(res)
load("low_dimension_original_result_1_40.rda")
processed_res_case2 <- pre_process_low(res)
load("low_dimension_original_result_1_60.rda")
processed_res_case3 <- pre_process_low(res)

res1 <- round(apply(processed_res_case1, c(1, 2), mean), 2)[c(1, 3), ]
res2 <- round(apply(processed_res_case2, c(1, 2), mean), 2)[c(1, 3), ]
res3 <- round(apply(processed_res_case3, c(1, 2), mean), 2)[c(1, 3), ]
rownames(res3) <- rownames(res2) <- rownames(res1) <- c("ABESS", "ASR-SIC")
colnames(res3) <- colnames(res2) <- colnames(res1) <- c("TPR", "TNR", "ReErr", "SLE")
res1
res2
res3

res1 <- round(apply(processed_res_case1, c(1, 2), sd), 2)[c(1, 3), ]
res2 <- round(apply(processed_res_case2, c(1, 2), sd), 2)[c(1, 3), ]
res3 <- round(apply(processed_res_case3, c(1, 2), sd), 2)[c(1, 3), ]
rownames(res3) <- rownames(res2) <- rownames(res1) <- c("ABESS", "ASR-SIC")
colnames(res3) <- colnames(res2) <- colnames(res1) <- c("TPR", "TNR", "ReErr", "SLE")
res1
res2
res3

######################################################################
############################### SI ###################################
######################################################################
pre_process <- function(res) {
  selected_metric <- c(5, 6, 8, 4, 1)
  selected_method <- 1:19
  res_tmp <- res[, selected_metric, ]
  res_tmp <- res_tmp[selected_method, , ]
  res_tmp[, 2, ] <- 1 - res_tmp[, 2, ]
  if (length(selected_metric) == 4) {
    res_tmp[, 3, ] <- res_tmp[, 3, ] - 10
  }
  res_tmp <- res_tmp[, c(1, 2, 4, 3, 5), ]
  res_tmp[, 3, ] <- res_tmp[, 3, ] * 10^3
  res_tmp
}

latex_max_index <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
latex_abs_index <- c(FALSE, FALSE, FALSE, TRUE, FALSE)
col_name <- c("Sensitivity", "Specificity", "ReErr", "SLE", "Runtimes")
row_name <- c("ABESS", "SDAR-EBIC", "SDAR-SIC", "CD-CV", "CD-SIC", "CDPSI-CV", "CDPSI-SIC", 
              "MCP-CV", "MCP-SIC", "MCP-EBIC", "MCP-BIC", 
              "SCAD-CV", "SCAD-SIC", "SCAD-EBIC", "SCAD-BIC", 
              "LASSO-CV", "LASSO-SIC", "LASSO-EBIC", "LASSO-BIC")

######################################################
################## None correlation ##################
######################################################
load("high_dimension_original_result_1_500.rda")
processed_res_p100 <- pre_process(res)
load("high_dimension_original_result_1_1500.rda")
processed_res_p500 <- pre_process(res)
load("high_dimension_original_result_1_2500.rda")
processed_res_p2500 <- pre_process(res)

sub_table1 <- apply(processed_res_p100, c(1, 2), mean)
rownames(sub_table1) <- NULL
sub_table2 <- apply(processed_res_p500, c(1, 2), mean)
rownames(sub_table2) <- NULL
sub_table3 <- apply(processed_res_p2500, c(1, 2), mean)
rownames(sub_table3) <- NULL
colnames(sub_table1) <- col_name
printbold(xtable(sub_table1, digits = 3), max = latex_max_index, abs = latex_abs_index)
printbold(xtable(sub_table2, digits = 3), max = latex_max_index, abs = latex_abs_index)
printbold(xtable(sub_table3, digits = 3), max = latex_max_index, abs = latex_abs_index)

sub_table1 <- apply(processed_res_p100, c(1, 2), sd)
sub_table2 <- apply(processed_res_p500, c(1, 2), sd)
sub_table3 <- apply(processed_res_p2500, c(1, 2), sd)
printbold(xtable(sub_table1, digits = 3), max = FALSE)
printbold(xtable(sub_table2, digits = 3), max = FALSE)
printbold(xtable(sub_table3, digits = 3), max = FALSE)

##########################################################
################## Constant correlation ##################
##########################################################
load("high_dimension_original_result_2_500.rda")
processed_res_p100 <- pre_process(res)
load("high_dimension_original_result_2_1500.rda")
processed_res_p500 <- pre_process(res)
load("high_dimension_original_result_2_2500.rda")
processed_res_p2500 <- pre_process(res)

sub_table1 <- apply(processed_res_p100, c(1, 2), mean)
rownames(sub_table1) <- NULL
sub_table2 <- apply(processed_res_p500, c(1, 2), mean)
rownames(sub_table2) <- NULL
sub_table3 <- apply(processed_res_p2500, c(1, 2), mean)
rownames(sub_table3) <- NULL
colnames(sub_table1) <- col_name
printbold(xtable(sub_table1, digits = 3), max = latex_max_index, abs = latex_abs_index)
printbold(xtable(sub_table2, digits = 3), max = latex_max_index, abs = latex_abs_index)
printbold(xtable(sub_table3, digits = 3), max = latex_max_index, abs = latex_abs_index)

sub_table1 <- apply(processed_res_p100, c(1, 2), sd)
sub_table2 <- apply(processed_res_p500, c(1, 2), sd)
sub_table3 <- apply(processed_res_p2500, c(1, 2), sd)
rownames(sub_table1) <- NULL
rownames(sub_table2) <- NULL
rownames(sub_table3) <- NULL
printbold(xtable(sub_table1, digits = 3), max = FALSE)
printbold(xtable(sub_table2, digits = 3), max = FALSE)
printbold(xtable(sub_table3, digits = 3), max = FALSE)

