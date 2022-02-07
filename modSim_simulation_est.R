# jonashaslbeck@gmail.com; March 2018

# Get iter from command line ...
iter <- commandArgs(trailingOnly=TRUE)

# -----------------------------------------------------------------------------
# ---------- Load Packages ----------------------------------------------------
# -----------------------------------------------------------------------------

library(devtools)

# Parallelize
library(foreach)
library(doParallel)

# Estimation
library(mgm) # and data generation
library(NetworkComparisonTest)
library(EstimateGroupNetwork)
library(gtools)
library(corpcor)


# -----------------------------------------------------------------------------
# ---------- Load Data ------------------------------------------------------
# -----------------------------------------------------------------------------


dataGen <- readRDS(file = paste0('output/ModSim_Data_Iter_', iter, '.RDS'))


# -----------------------------------------------------------------------------
# ---------- Estimation -------------------------------------------------------
# -----------------------------------------------------------------------------

# --------- Global Settings ---------

n_seq <- round(exp(seq(3.4, 7.5, length=12))) # n-variations
p <- 13 # number of variables
gamma <- .5 # EBIC tuning parameter


# --------- Parallelize estimation over n-variations ---------

nClust <- 12 # We run the n-variations in parallel

# Fire up foreach
cl <- makeCluster(nClust)
registerDoParallel(cl)


# start clustered loops
outlist <- foreach(n=1:12,
                   .packages=c("mgm", "NetworkComparisonTest", "EstimateGroupNetwork"),
                   .verbose=FALSE) %dopar% {
                     
                     
                     # ----- Subset Data -----
                     
                     data <- dataGen$Datalist$data_finite[1:n_seq[n], ]
                     type <- rep("g", p)
                     level <- rep(1, p)
                     
                     
                     # ----- Create Output List -----
                     
                     est_list <- list("MNW_1_OR" = NULL,
                                      "MNW_2_OR" = NULL,
                                      "MNW_3_OR" = NULL,
                                      "MNW_1_AND" = NULL,
                                      "MNW_2_AND" = NULL,
                                      "MNW_3_AND" = NULL,
                                      "NCT_1" = NULL,
                                      "NCT_2" = NULL,
                                      "FGL_1" = NULL,
                                      "FGL_2" = NULL,
                                      "Runtime" = NULL)
                     
                     # start timer
                     time_start <- proc.time()[3]
                     
                     # -------------------
                     # ----- MODELS ------
                     # -------------------
                     
                     # ----- Moderated Network Models (1); moderators = m ------                     
                     
                     # Model Selection with EBIC + OR-rule
                     set.seed(iter)
                     MNW_1_OR <- mgm(data = data,
                                     type = type,
                                     level = level,
                                     k = 2,
                                     moderators = dataGen$m,
                                     lambdaSel = "EBIC",
                                     lambdaGam = gamma,
                                     ruleReg = "OR", 
                                     pbar = FALSE)
                     
                     est_list$MNW_1_OR <- MNW_1_OR
                     
                     
                     # Model Selection with EBIC + AND-rule                     
                     set.seed(iter)
                     MNW_1_AND <- mgm(data = data,
                                      type = type,
                                      level = level,
                                      k = 2,
                                      moderators = dataGen$m,
                                      lambdaSel = "EBIC",
                                      lambdaGam = gamma,
                                      ruleReg = "AND", 
                                      pbar = FALSE)
                     
                     est_list$MNW_1_AND <- MNW_1_AND
                     
                     
                     # ----- Moderated Network Models (2); moderators = 1:p, sequentially ------                     
                     
                     
                     # Model Selection with EBIC + OR-rule
                     MNW_2_OR <- list()
                     
                     for(i in 1:p) {
                       
                       set.seed(iter)
                       MNW_2_OR[[i]] <- mgm(data = data,
                                            type = type,
                                            level = level,
                                            k = 2,
                                            moderators = i,
                                            lambdaSel = "EBIC",
                                            lambdaGam = gamma,
                                            ruleReg = "OR", 
                                            pbar = FALSE)
                       
                     } # end for i (moderators)
                     
                     est_list$MNW_2_OR <- MNW_2_OR
                     
                     
                     # Model Selection with EBIC + AND-rule      
                     
                     MNW_2_AND <- list()
                     
                     for(i in 1:p) {
                       MNW_2_AND[[i]] <- mgm(data = data,
                                             type = type,
                                             level = level,
                                             k = 2,
                                             moderators = i,
                                             lambdaSel = "EBIC",
                                             lambdaGam = gamma,
                                             ruleReg = "AND", 
                                             pbar = FALSE)
                     } # end for i (moderators)
                     
                     est_list$MNW_2_AND <- MNW_2_AND
                     
                     
                     # ----- Moderated Network Models (3); k = 3 ------                     
                     
                     # Model Selection with EBIC + OR-rule
                     
                     MNW_3_OR <- mgm(data = data,
                                     type = type,
                                     level = level,
                                     k = 3,
                                     lambdaSel = "EBIC",
                                     lambdaGam = gamma,
                                     ruleReg = "OR", 
                                     pbar = FALSE)
                     
                     est_list$MNW_3_OR <- MNW_3_OR
                     
                     
                     # Model Selection with EBIC + AND-rule                           
                     
                     MNW_3_AND <- mgm(data = data,
                                      type = type,
                                      level = level,
                                      k = 3,
                                      lambdaSel = "EBIC",
                                      lambdaGam = gamma,
                                      ruleReg = "AND", 
                                      pbar = FALSE)
                     
                     est_list$MNW_3_AND <- MNW_3_AND
                     
                     
                     # ----- NCT (1); moderators = m ------                     
                     
                     # Split half
                     data_ord <- data[order(data[, dataGen$m]), ]
                     md_data <- median(data[, dataGen$m])
                     
                     data_0 <- data_ord[data_ord[, dataGen$m] < md_data, ]
                     data_1 <- data_ord[data_ord[, dataGen$m] >= md_data, ]
                     
                     # Run NCT
                     NCT_1 <- NCT(data1 = data_0, 
                                  data2 = data_1, 
                                  gamma = .5, # Default value of NCT
                                  it = 1000, # run 1000 iterations
                                  test.edges = TRUE, 
                                  binary.data = FALSE,
                                  edges = "all",
                                  progressbar = FALSE)
                     
                     est_list$NCT_1 <- NCT_1
                     
                     
                     
                     # ----- NCT (2); moderators = 1:p, sequentially ------                     
                     
                     
                     NCT_2 <- list()
                     
                     for(i in 1:p) {
                       
                       set.seed(1)
                       
                       # Split half
                       data_ord <- data[order(data[, i]), ]
                       md_data <- median(data[, i])
                       
                       data_0 <- data_ord[data_ord[, i] < md_data, ]
                       data_1 <- data_ord[data_ord[, i] >= md_data, ]
                       
                       # Fit NCT
                       NCT_2[[i]] <- NCT(data1 = data_0,
                                         data2 = data_1,
                                         gamma = .5, # Default value of NCT
                                         it = 1000, # run 1000 iterations
                                         test.edges = TRUE,
                                         binary.data = FALSE,
                                         edges = "all",
                                         progressbar = FALSE)
                       
                     }
                     
                     est_list$NCT_2 <- NCT_2
                     
                     # ----- FGL (1); moderators = m ------    
                     
                     # Split half
                     data_ord <- data[order(data[, dataGen$m]), ]
                     md_data <- median(data[, dataGen$m])
                     
                     data_0 <- data_ord[data_ord[, dataGen$m] < md_data, ]
                     data_1 <- data_ord[data_ord[, dataGen$m] >= md_data, ]
                     
                     l_data <- list("part0" = data_0, 
                                    "part1" = data_1)
                     
                     # Fit FGL
                     FGL_1 <- EstimateGroupNetwork(X = l_data, 
                                                   inputType = "list.of.dataframes", 
                                                   method = "InformationCriterion", 
                                                   strategy = "sequential",
                                                   gamma = .5)
                     
                     
                     est_list$FGL_1 <- FGL_1
                     
                     
                     # ----- FGL (2); moderators = 1:p, sequentially ------
                     
                     
                     FGL_2 <- list()
                     
                     for(i in 1:p) {
                       
                       # Split half
                       data_ord <- data[order(data[, i]), ]
                       md_data <- median(data[, i])
                       
                       data_0 <- data_ord[data_ord[, i] < md_data, ]
                       data_1 <- data_ord[data_ord[, i] >= md_data, ]
                       
                       l_data <- list("part0" = data_0, 
                                      "part1" = data_1)
                       
                       
                       # Fit FGL
                       FGL_2[[i]] <- EstimateGroupNetwork(X = l_data, 
                                                          inputType = "list.of.dataframes", 
                                                          method = "InformationCriterion", 
                                                          strategy = "sequential",
                                                          gamma = .5)
                       
                       
                     }
                     
                     est_list$FGL_2 <- FGL_2
                     
                     
                     # ----- Return -----
                     
                     est_list$Runtime <- proc.time()[3] - time_start
                     
                     return(est_list)
                     
                   } # end: foreach (n)


# Close down foreach
stopCluster(cl)


# --------- Save Results ---------

saveRDS(outlist, file = paste0('ModSim_Est_Iter_', iter, '.RDS'))



