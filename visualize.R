rm(list = ls()); gc(reset = TRUE)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(abind)

result_path <- "./simulation_result"
setwd(result_path)

####################################################################
######################### Helpful functions ########################
####################################################################
## For simulation result preprocessing:
pre_process_low <- function(res) {
  selected_metric <- c(5, 6, 8, 4)
  selected_method <- 1:3
  res_tmp <- res[, selected_metric, ]
  res_tmp <- res_tmp[selected_method, , ]
  res_tmp[, 2, ] <- 1 - res_tmp[, 2, ]
  res_tmp[, 3, ] <- (res_tmp[, 3, ] - 3)
  res_tmp
}

pre_process <- function(res) {
  selected_metric <- c(5, 6, 8, 4)
  selected_method <- c(1, 8, 9, 12, 13, 16, 17)
  res_tmp <- res[, selected_metric, ]
  res_tmp <- res_tmp[selected_method, , ]
  res_tmp[, 2, ] <- 1 - res_tmp[, 2, ]
  if (length(selected_metric) == 4) {
    res_tmp[, 3, ] <- res_tmp[, 3, ] - 10
  }
  res_tmp
}

runtime_pre_process <- function(res) {
  res_tmp <- res
  selected_metric <- c(1)
  res_tmp <- res_tmp[, selected_metric, , drop = TRUE]
  selected_method <- c(1, 2)
  res_tmp <- res_tmp[selected_method, ]
  rownames(res_tmp) <- NULL
  res_tmp <- as.data.frame(res_tmp)
  res_tmp[["method"]] <- factor(c("ABESS", "ASR"), 
                                levels = c("ABESS", "ASR"))
  
  res_tmp
}

performance_pre_process <- function(res) {
  res_tmp <- res
  selected_metric <- c(5, 6, 4)
  res_tmp <- res_tmp[, selected_metric, , drop = TRUE]
  res_tmp[, 2, ] <- 1 - res_tmp[, 2, ]
  selected_method <- c(1, 3)
  res_tmp <- res_tmp[selected_method, , ]
  res_tmp <- apply(res_tmp, MARGIN = c(1, 2), mean)
  rownames(res_tmp) <- NULL
  colnames(res_tmp) <- c("TPR", "TNR", "ReErr")
  res_tmp <- as.data.frame(res_tmp)
  res_tmp[1, ] <- res_tmp[1, ] - res_tmp[2, ]
  res_tmp <- res_tmp[1, ]
  res_tmp <- melt(res_tmp)
  res_tmp
}

runtime_pre_process_high <- function(res) {
  res_tmp <- res
  selected_metric <- c(1)
  res_tmp <- res_tmp[, selected_metric, , drop = TRUE]
  selected_method <- c(1, 8, 9, 12, 13, 16, 17)
  res_tmp <- res_tmp[selected_method, ]
  rownames(res_tmp) <- NULL
  res_tmp <- as.data.frame(res_tmp)
  row_name <- c("ABESS", "MCP-CV", "MCP-SIC",
                "SCAD-CV", "SCAD-SIC", "LASSO-CV", "LASSO-SIC")

  res_tmp[["method"]] <- factor(row_name, levels = row_name)
  res_tmp
}

## For Figure display:
figure_display_batch <- function(res) {
  res_tmp <- res
  size <- dim(res_tmp)[3]
  res_list <- list()
  for (i in 1:size) {
    res_list[[i]] <- res_tmp[, , i]
  }
  
  res_list <- do.call("rbind", res_list)
  if (ncol(res_list) == 4) {
    col_name <- c("TPR", "TNR", "SLE", "ReErr")
  } else {
    col_name <- c("TPR", "TNR", "ReErr")
  }
  colnames(res_list) <- col_name
  res_list <- as.data.frame(res_list)
  rownames(res_list) <- NULL
  row_name <- c("ABESS", "MCP-CV", "MCP-SIC", "SCAD-CV", "SCAD-SIC", "LASSO-CV", "LASSO-SIC")
  res_list[["method"]] <- rep(row_name, size)
  res_list[["method"]] <- factor(res_list[["method"]], 
                                 levels = row_name)
  res_list[["p"]] <- rep(c(500, 1500, 2500), each = length(row_name) * size / 3)
  res_list[["p"]] <- as.factor(res_list[["p"]])
  p_dat <- melt(res_list, id.vars = c("p", "method"))
  factor_level <- c("TPR", "TNR", "ReErr", "SLE")
  p_dat[["variable"]] <- factor(p_dat[["variable"]], 
                                levels = factor_level)
  
  p <- ggplot(p_dat, aes(x = p, y = value, fill = method)) + 
    geom_boxplot() +
    facet_wrap(. ~ variable, scales = "free", ncol = 1) + 
    xlab("dimension") + 
    ylab("") + 
    theme_bw() + 
    theme(legend.position = "bottom", 
          legend.margin = margin(-15, -20, -3, -20),
          legend.box.margin = margin(15, 15, 15, 15)) + 
    guides(fill = guide_legend(byrow = FALSE))
  p
}

performance_figure_display_batch <- function(res) {
  p_dat <- res
  factor_level <- c("TPR", "TNR", "ReErr")
  p_dat[["variable"]] <- factor(p_dat[["variable"]], 
                                levels = factor_level)
  p <- ggplot(p_dat, aes(x = dimensionality, y = value, colour = variable)) +
    geom_point() + geom_line() + theme_bw() + 
    scale_y_continuous(breaks = (-5:5) * 1e-3, limits = c(-0.005, 0.005)) +
    scale_color_discrete("Metric") + 
    theme(legend.position = c(0.5, 0.9), 
          legend.box = "horizontal", 
          legend.direction = "horizontal", 
          legend.box.background = element_rect(colour = "black")) + 
    ylab("difference") + xlab("dimension")
  p
}

runtime_figure_display_low <- function(res, formula, algorithm_type = "ABESS", 
                                       color = "#2E9FDF", xy_coord = c(1, 2), decimal = 4) {
  p_dat <- melt(res, id.vars = c("dimensionality"), value.name = "runtime")
  
  if (decimal == 4) {
    scaleFUN <- function(x) sprintf("%.4f", x)
  } else {
    scaleFUN <- function(x) sprintf("%.3f", x)
  }
  
  y_axis_name <- paste0(algorithm_type, "'s runtime (sec)")
  x_axis_name <- "dimension"
  
  p <- ggplot(p_dat, aes(x = dimensionality, y = runtime)) +
    geom_point(color = color, size = 0.5) +
    geom_jitter(width = 0.15, color = color, size = 0.5) +
    geom_smooth(colour = color, formula = formula, method = "lm", se = FALSE) + 
    theme(legend.position = "none") + 
    scale_y_continuous(labels = scaleFUN) + 
    ylab(y_axis_name) + xlab(x_axis_name) + 
    theme_bw()
  p
}

runtime_figure_display_high <- function(res) {
  p_dat <- melt(res, id.vars = c("dimensionality", "method", "type"))
  
  p_dat <- p_dat %>% 
    group_by(type, dimensionality, method) %>% 
    summarise(runtime = mean(value), runtime_sd = sd(value), 
              runtime_minus_1se = runtime - runtime_sd, 
              runtime_add_1se = runtime + runtime_sd) 
  
  pd <- position_dodge(56)
  
  p <- ggplot(p_dat, 
              aes(x = dimensionality, y = runtime, 
                  ymin = runtime_minus_1se, ymax = runtime_add_1se, 
                  color = method)) + 
    facet_wrap(. ~ type) + 
    geom_errorbar(position = pd, width = 188) +
    geom_point(position = pd) + 
    geom_line(position = pd) + 
    theme_bw() + 
    ylab("runtime (sec)") + 
    xlab("dimension") + 
    theme(legend.position = "bottom", 
          legend.margin = margin(-18, -20, -3, -20),
          legend.box.margin = margin(15, 15, 15, 15))
  p
}

#########################################################################
######################### Figure 1 (low dimension) ######################
#########################################################################
for (p in 20:40) {
  load(paste0("runtime_original_result_", p, ".rda"))
  dimensionality <- p
  if (p == 20) {
    runtime_processed_res <- cbind.data.frame(runtime_pre_process(res), "dimensionality" = dimensionality)
    performance_processed_res <- cbind.data.frame(performance_pre_process(res), "dimensionality" = dimensionality)
  } else {
    runtime_processed_res1 <- cbind.data.frame(runtime_pre_process(res), "dimensionality" = dimensionality)
    runtime_processed_res <- rbind(runtime_processed_res, runtime_processed_res1)
    performance_processed_res1 <- cbind.data.frame(performance_pre_process(res), "dimensionality" = dimensionality)
    performance_processed_res <- rbind(performance_processed_res, performance_processed_res1)
  }
}

p1 <- performance_figure_display_batch(performance_processed_res)
p1

mycols <- c("#2E9FDF", "#FC4E07")
p2_dat <- runtime_processed_res[runtime_processed_res[["method"]] == "ABESS", ]
p2_dat[["method"]] <- NULL
formula <- y ~ x
p2 <- runtime_figure_display_low(p2_dat, formula, "ABESS", mycols[1], xy_coord = c(27, 6.5*1e-3))
p2

p3_dat <- runtime_processed_res[runtime_processed_res[["method"]] == "ASR", ]
p3_dat[["method"]] <- NULL
formula <- y ~ I(2^x)
p3 <- runtime_figure_display_low(p3_dat, formula, "ASR", mycols[2], c(27, 15), decimal = 3)
p3

p <- ggarrange(p1, p2, p3, labels = c("A", "B", "C"), ncol = 1, 
               font.label = list(size = 14, color = "black", 
                                 face = "plain", family = NULL))
p
ggexport(p, filename = "runtime.pdf", height = 9, width = 6)

###################################################################
##################### Figure 2 (high dimension) ###################
###################################################################
load("high_dimension_original_result_1_500.rda")
processed_res_p100 <- pre_process(res)
load("high_dimension_original_result_1_1500.rda")
processed_res_p500 <- pre_process(res)
load("high_dimension_original_result_1_2500.rda")
processed_res_p2500 <- pre_process(res)
processed_res <- abind(processed_res_p100, processed_res_p500, processed_res_p2500)
p1 <- figure_display_batch(processed_res)

load("high_dimension_original_result_2_500.rda")
processed_res_p100 <- pre_process(res)
load("high_dimension_original_result_2_1500.rda")
processed_res_p500 <- pre_process(res)
load("high_dimension_original_result_2_2500.rda")
processed_res_p2500 <- pre_process(res)
processed_res <- abind(processed_res_p100, processed_res_p500, processed_res_p2500)
p2 <- figure_display_batch(processed_res)

p <- ggarrange(p1, p2, labels = c("A", "B"), ncol = 2, hjust = c(-2, -2), 
               common.legend = TRUE, legend = "bottom", 
               font.label = list(size = 14, color = "black", 
                                 face = "plain", family = NULL))
p
ggsave(p, filename = "high_dimension_result_1.pdf", height = 9.8, width = 7.6)

for (type in c(1, 2)) {
  for (p in c(1, 3, 5) * 500) {
    load(paste0("high_dimension_original_result_", type, "_", p, ".rda"))
    type_str <- ifelse(type == 1, "Independence", "Constant")
    if (p == 500 & type == 1) {
      runtime_all <- cbind.data.frame(runtime_pre_process_high(res), "dimensionality" = p, "type" = type_str)
    } else {
      runtime_tmp <- cbind.data.frame(runtime_pre_process_high(res), "dimensionality" = p, "type" = type_str)
      runtime_all <- rbind(runtime_all, runtime_tmp)
    }
  }
}

p3 <- runtime_figure_display_high(runtime_all)
p3

p <- ggarrange(p3, labels = c("C"), ncol = 1, 
               font.label = list(size = 14, color = "black", 
                                 face = "plain", family = NULL))
p
ggsave(p, filename = "high_dimension_result_2.pdf", height = 3.8, width = 6.8)


setwd("..")