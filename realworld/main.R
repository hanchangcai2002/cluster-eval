rm(list = ls())
packages <- c("dplyr", "lcmm", "mixAK","gbmt", "flexmix", "mclust", "clusterCrit")
to_install <- setdiff(packages, rownames(installed.packages()))
print(paste("Installing packages:", to_install))
if (length(to_install)) install.packages(to_install)
invisible(lapply(packages, library, character.only = TRUE))


setwd("/ocean/projects/med220007p/hcai5/jinyuan/out/realworld")
source("/ocean/projects/med220007p/hcai5/jinyuan/code/realworld_func_0920.R")


realworld_res <- list()
dat <- readRDS("./realworld.RDS")
realworld_res[["realworld_flexmix"]] <- run_flexmix_analysis(dat, "realworld", c("sbp", "dbp"))
realworld_res[["realworld_gbmt"]] <- run_gbmt_analysis(dat, "realworld", c("sbp", "dbp"))
realworld_res[["realworld_lcmm"]] <- run_lcmm_analysis(dat, "realworld", c("sbp", "dbp"))
realworld_res[["realworld_mixAK"]] <- run_mixAK_analysis(dat, "realworld", c("sbp", "dbp"))

saveRDS(realworld_res, "./Sep22_realworld_res.RDS")