# jonashaslbeck@gmail.com; March 2018

# Get iter from command line ...
iter <- commandArgs(trailingOnly=TRUE)

# -----------------------------------------------------------------------------
# ---------- Load Packages ----------------------------------------------------
# -----------------------------------------------------------------------------

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
# ---------- Generate Graph ---------------------------------------------------
# -----------------------------------------------------------------------------

dataGen <- list("Graph" = NULL,
                "Datalist" = NULL)

# A) ----------- Specify Graph ----------

# ----- Function to sample 1 graph -----

GenerateGraph <- function() {
  p <- 13 # number of variables 
  p_e <- p-1
  n_edges <- (p_e*(p_e-1)/2)
  select_edges <- sample(1:n_edges, size = 6, replace = FALSE) # obtain edges
  select_types <- sample(c(1, 1, 2, 2, 3, 3), size = 6, replace = FALSE) # obtain edge-types
  v_lower_tri <- rep(0, n_edges)
  v_lower_tri[select_edges] <- select_types
  
  # Put all in graph
  G <- matrix(0, p, p)
  Gsm <- matrix(0, p_e, p_e)
  Gsm[lower.tri(Gsm)] <- v_lower_tri
  
  # Constraint I: no (additional) edges to moderator
  G[1:p_e, 1:p_e] <- Gsm # node 13 is always the moderator
  G <- G + t(G)
  return(G)
}


# ----- Sample graphs until constraints on graph structure are satisfied -----

set.seed(iter)

check <- 2
counter <- 0

while(check > 0) {
  
  G <- GenerateGraph()
  
  # Constraint II: Degree < 3 for all nodes
  G_bn <- G
  G_bn[G_bn != 0] <- 1
  if(any(colSums(G_bn) > 2)) check <- 1 else check <- 0
  
  # Avoid invite loop
  if(counter > 1000) stop("Graph resampling failed!")
  if(check == 0) break
  counter <- counter + 1
  
}

dataGen$graphResample <- counter # save resampling steps

# Obtain Moderator variable
dataGen$m <- m <- 13 # fix moderator

# get interactions PW
G_mod <- G
G_mod[upper.tri(G_mod)] <- 0 # to avoid getting each interaction twice
int_type_1 <- which(G_mod == 1, arr.ind = TRUE)
int_type_2 <- which(G_mod == 2, arr.ind = TRUE)
int_type_3 <- which(G_mod == 3, arr.ind = TRUE)


# B) ----------- Input for mgmsampler() -----------

# a) General Graph Info
p <- 13
type = rep('g', p) 
level = rep(1, p) 
theta <- .2

# b) Define interactions
factors <- list()
interactions <- list()

# b.1) Pairwise Effects (4)
factors[[1]] <- matrix(c(int_type_1[1, ],
                         int_type_1[2, ],
                         int_type_2[1, ],
                         int_type_2[2, ]), ncol=2, byrow = T) # 4 pairwise interactions (type 1 and 2)

interactions[[1]] <- vector('list', length = 3)
interactions[[1]][[1]] <- array(theta, dim=c(1, 1))
interactions[[1]][[2]] <- array(theta, dim=c(1, 1))
interactions[[1]][[3]] <- array(theta, dim=c(1, 1))
interactions[[1]][[4]] <- array(theta, dim=c(1, 1))


# b.2) Moderation Effects (4)
factors[[2]] <- matrix(c(int_type_2[1, ], m,
                         int_type_2[2, ], m,
                         int_type_3[1, ], m,
                         int_type_3[2, ], m), ncol=3, byrow = T) # 4 pairwise interactions (type 2 and 3)
interactions[[2]] <- vector('list', length = 3)
interactions[[2]][[1]] <- array(theta, dim=c(1, 1, 1))
interactions[[2]][[2]] <- array(theta, dim=c(1, 1, 1))
interactions[[2]][[3]] <- array(theta, dim=c(1, 1, 1))
interactions[[2]][[4]] <- array(theta, dim=c(1, 1, 1))


# c) Define Thresholds
thresholds <- vector('list', length = p)
thresholds <- lapply(thresholds, function(x) 0 ) # all intercepts are zero

# d) Define Variances
sds <- rep(1, p) # All variances equal to 1

dataGen$Graph <- G


# -----------------------------------------------------------------------------
# ---------- Sample Data ------------------------------------------------------
# -----------------------------------------------------------------------------

set.seed(iter)

time_start <- proc.time()[3]
dlist <- mgmsampler(factors = factors,
                    interactions = interactions,
                    thresholds = thresholds,
                    sds = sds,
                    type = type,
                    level = level,
                    N = 10000,
                    nIter = 100,
                    pbar = FALSE)
time_end <- proc.time()[3] - time_start

# Subset cases in 99.9% quantile
data_finite <- dlist$data
data_finite <- data_finite[apply(data_finite, 1, function(x) all(abs(x) < qnorm(c(.999), 0, 1))),]
data_finite <- data_finite[apply(data_finite, 1, function(x) all(!is.na(x))), ]

# Save sampling info in data list
dlist$pop_finite <- nrow(data_finite) / nrow(dlist$data) # save: how many finite?
dlist$data_finite <- data_finite
dlist$time_sampling <- time_end

dataGen$Datalist <- dlist

# Save Data & Graph
saveRDS(dataGen, file = paste0('ModSim_Data_Iter_', iter, '.RDS'))

if(nrow(data_finite) < 1808) stop("Convergence rate too low for estimation!")




