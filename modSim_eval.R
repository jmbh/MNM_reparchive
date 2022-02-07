# jonashaslbeck@gmail.com; February 2018 

# -----------------------------------------------------------------
# ----------------- Load Results ----------------------------------
# -----------------------------------------------------------------


# !!! Specity Directories of Simulation Output !!!

outdir_data <- ".../output_data/" # Specify directory of data output of simulation
outdir_est <- ".../output_est/" # Specify directory of estimation output of simulation


first100 <- readRDS(file = "data/first100.RDS") # Subset 100 cases with largest proportion of finite values (created in ChecksOnData.R)

# Load Data files
dfiles_data <- list.files(outdir_data)
iters <- substr(dfiles_data, 18, 20) # get iters
iters <- sub(".", "", iters, fixed=TRUE)
iters <- sub("R", "", iters, fixed=TRUE)
iters <- as.numeric(iters)
ind_in_first100 <- iters %in% first100
dfiles_data <- dfiles_data[ind_in_first100]

l_data <- list()
for(i in 1:length(dfiles_data)) {
  l_data[[i]] <- readRDS(paste0(outdir_data, dfiles_data[i]))
  print(i)
} 

# Load Estimation files
dfiles_est <- list.files(outdir_est)
iters <- substr(dfiles_est, 17, 19) # get iters
iters <- sub(".", "", iters, fixed=TRUE)
iters <- sub("R", "", iters, fixed=TRUE)
iters <- as.numeric(iters)
ind_in_first100 <- iters %in% first100
dfiles_est <- dfiles_est[ind_in_first100]

iters[ind_in_first100]

l_est <- list()
for(i in 1:length(dfiles_est)) {
  l_est[[i]] <- readRDS(paste0(outdir_est, dfiles_est[i]))
  print(i)
} 


# -----------------------------------------------------------------
# ----------------- Aux Functions ---------------------------------
# -----------------------------------------------------------------

source("f_EvalModNWs.R")

# -----------------------------------------------------------------
# ----------------- Preprocess ------------------------------------
# -----------------------------------------------------------------

# Output:
# - for Pairwise interactions: sensitivity/precision of the 3 ModNW algorithms
# - for 3-way interactions: sensitivity/precision of the 3 ModNW algorithms and the two NCT/FGL algorithms

nIter <- length(dfiles_data)
p <- 13
## Storage
# iteraction types, types, estimators, n-variation, Pfmeasure (Sen/Pre); 
#     interaction types: PW, PW with mod, mod with PW, mod 
a_res <- array(NA, dim = c(nIter, 4, 7, 12, 2)) 

ORRULE <- FALSE

for(iter in 1:nIter) {
  for(n in 1:12) {
    
    print(paste("iter", iter, "n", n))
    
    # ----- Find the interaction types -----
    
    G_mod <- l_data[[iter]]$Graph
    G_mod[upper.tri(G_mod)] <- 0 # to avoid getting each interaction twice
    
    int_type_1 <- which(G_mod == 1, arr.ind = TRUE)
    int_type_2 <- which(G_mod == 2, arr.ind = TRUE)
    int_type_3 <- which(G_mod == 3, arr.ind = TRUE)
    
    # Order to make comparisons below easier
    int_type_1 <- t(apply(int_type_1, 1, sort))
    int_type_2 <- t(apply(int_type_2, 1, sort))
    int_type_3 <- t(apply(int_type_3, 1, sort))
    
    m <- l_data[[iter]]$m # the moderator
    
    
    # ---------------------------------------------------
    # ----- Different Extraction for each Estimator ----- 
    # ---------------------------------------------------
    
    # ----- Moderated Network Models (1); moderators = m ------
    
    if(ORRULE) {
      est_2way <- l_est[[iter]][[n]]$MNW_1_OR$interactions$indicator[[1]]
      est_3way <- l_est[[iter]][[n]]$MNW_1_OR$interactions$indicator[[2]]  
    } else {
      est_2way <- l_est[[iter]][[n]]$MNW_1_AND$interactions$indicator[[1]]
      est_3way <- l_est[[iter]][[n]]$MNW_1_AND$interactions$indicator[[2]]  
    }
    
    
    outlist <- f_EvalModNWs(est_2way = est_2way,
                            est_3way = est_3way, 
                            int_type_1 = int_type_1, 
                            int_type_2 = int_type_2, 
                            int_type_3 = int_type_3, 
                            m = m)
    
    # Save in array
    a_res[iter, 1, 1, n, 1] <- outlist$type1_S
    a_res[iter, 2, 1, n, 1] <- outlist$type2_S
    a_res[iter, 1, 1, n, 2] <- outlist$type1_P
    a_res[iter, 2, 1, n, 2] <- outlist$type2_P
    
    a_res[iter, 3, 1, n, 1] <- outlist$type3_S
    a_res[iter, 4, 1, n, 1] <- outlist$type4_S
    a_res[iter, 3, 1, n, 2] <- outlist$type3_P
    a_res[iter, 4, 1, n, 2] <- outlist$type4_P
    
    
    # ----- Moderated Network Models (2); moderators = 1:p, sequentially ------                     
    
    # Collapse Estimates across p=10 runs
    est_2way <- est_3way <- list()
    
    if(ORRULE) {
      for(i in 1:10) {
        est_2way[[i]] <- l_est[[iter]][[n]]$MNW_2_OR[[i]]$interactions$indicator[[1]]
        est_3way[[i]] <- l_est[[iter]][[n]]$MNW_2_OR[[i]]$interactions$indicator[[2]]
      }
    } else {
      for(i in 1:10) {
        est_2way[[i]] <- l_est[[iter]][[n]]$MNW_2_AND[[i]]$interactions$indicator[[1]]
        est_3way[[i]] <- l_est[[iter]][[n]]$MNW_2_AND[[i]]$interactions$indicator[[2]]
      }
    }
    
    
    est_2way <- do.call(rbind, est_2way)
    est_3way <- do.call(rbind, est_3way)
    
    est_2way_sort <- t(apply(est_2way, 1, sort))
    est_3way_sort <- t(apply(est_3way, 1, sort))
    
    est_2way_final <- est_2way_sort[!duplicated(est_2way_sort), ]
    est_3way_final <- est_3way_sort[!duplicated(est_3way_sort), ]
    
    outlist <- f_EvalModNWs(est_2way = est_2way_final,
                            est_3way = est_3way_final, 
                            int_type_1 = int_type_1, 
                            int_type_2 = int_type_2, 
                            int_type_3 = int_type_3, 
                            m = m)
    
    # Save in array
    a_res[iter, 1, 2, n, 1] <- outlist$type1_S
    a_res[iter, 2, 2, n, 1] <- outlist$type2_S
    a_res[iter, 1, 2, n, 2] <- outlist$type1_P
    a_res[iter, 2, 2, n, 2] <- outlist$type2_P
    a_res[iter, 3, 2, n, 1] <- outlist$type3_S
    a_res[iter, 4, 2, n, 1] <- outlist$type4_S
    a_res[iter, 3, 2, n, 2] <- outlist$type3_P
    a_res[iter, 4, 2, n, 2] <- outlist$type4_P
    
    
    # ----- Moderated Network Models (3); k = 3 ------  
    
    if(ORRULE) {
      est_2way <- l_est[[iter]][[n]]$MNW_3_OR$interactions$indicator[[1]]
      est_3way <- l_est[[iter]][[n]]$MNW_3_OR$interactions$indicator[[2]]  
    } else {
      est_2way <- l_est[[iter]][[n]]$MNW_3_AND$interactions$indicator[[1]]
      est_3way <- l_est[[iter]][[n]]$MNW_3_AND$interactions$indicator[[2]]  
    }
    
    outlist <- f_EvalModNWs(est_2way = est_2way,
                            est_3way = est_3way, 
                            int_type_1 = int_type_1, 
                            int_type_2 = int_type_2, 
                            int_type_3 = int_type_3, 
                            m = m)
    
    # Save in array
    a_res[iter, 1, 3, n, 1] <- outlist$type1_S
    a_res[iter, 2, 3, n, 1] <- outlist$type2_S
    a_res[iter, 1, 3, n, 2] <- outlist$type1_P
    a_res[iter, 2, 3, n, 2] <- outlist$type2_P
    a_res[iter, 3, 3, n, 1] <- outlist$type3_S
    a_res[iter, 4, 3, n, 1] <- outlist$type4_S
    a_res[iter, 3, 3, n, 2] <- outlist$type3_P
    a_res[iter, 4, 3, n, 2] <- outlist$type4_P
    
    
    # ----- NCT (1); moderators = m ------           
    
    M <- l_est[[iter]][[n]]$NCT_1 # NCT object
    
    pw_sig <- as.matrix(M$einv.pvals[M$einv.pvals[,3] < .05, ])
    n_sig <- nrow(pw_sig)
    M_comp <- matrix(NA, nrow = n_sig, ncol = 8)
    M_comp[, 1:3] <- pw_sig
    
    if(nrow(pw_sig) == 0) {
      
      a_res[iter, 3, 4, n, 1] <- 0
      a_res[iter, 4, 4, n, 1] <- 0
      a_res[iter, 3, 4, n, 2] <- NA
      a_res[iter, 4, 4, n, 2] <- NA
      
      
    } else {
      
      # get sign of difference
      for(s in 1:n_sig) M_comp[s, 4] <- sign(M$nw2[matrix(pw_sig[s,1:2], ncol=2)] -  M$nw1[matrix(pw_sig[s,1:2], ncol=2)])
      
      # Sensitivity for Mod (with Pw)
      M_comp[, 5:6] <- matrix(rep(int_type_2[1, ], n_sig), ncol=2, byrow = TRUE)
      M_comp[, 7:8] <- matrix(rep(int_type_2[2, ], n_sig), ncol=2, byrow = TRUE)
      v_pwMod_sen <- apply(M_comp, 1, function(x) x[4] == 1 & (identical(x[1:2], x[5:6]) | identical(x[1:2], x[7:8])) )
      
      # Sensitivty for Mod (alone)
      M_comp[, 5:6] <- matrix(rep(int_type_3[1, ], n_sig), ncol=2, byrow = TRUE)
      M_comp[, 7:8] <- matrix(rep(int_type_3[2, ], n_sig), ncol=2, byrow = TRUE)
      v_Mod_sen <- apply(M_comp, 1, function(x) x[4] == 1 & (identical(x[1:2], x[5:6]) | identical(x[1:2], x[7:8])) )
      
      
      a_res[iter, 3, 4, n, 1] <- sum(v_pwMod_sen) / 2 
      a_res[iter, 4, 4, n, 1] <- sum(v_Mod_sen) / 2
      a_res[iter, 3, 4, n, 2] <- (sum(v_pwMod_sen) + sum(v_Mod_sen)) / n_sig
      a_res[iter, 4, 4, n, 2] <- (sum(v_pwMod_sen) + sum(v_Mod_sen)) / n_sig
      
    }
    
    
    
    
    # ----- NCT (2); moderators = 1:p, sequentially ------  
    
    # Combine Results
    
    l_M <- l_est[[iter]][[n]]$NCT_2 # NCT objects
    l_M_comp <- list()
    
    for(i in 1:p) {
      
      M <- l_M[[i]]
      
      pw_sig <- as.matrix(M$einv.pvals[M$einv.pvals[,3] < .05, ])
      n_sig <- nrow(pw_sig)
      
      M_comp <- matrix(NA, nrow = n_sig, ncol = 8)
      
      if(n_sig != 0) {
        
        M_comp[, 1:3] <- pw_sig
        # get sign of difference
        for(s in 1:n_sig) M_comp[s, 4] <- sign(M$nw2[matrix(pw_sig[s,1:2], ncol=2)] -  M$nw1[matrix(pw_sig[s,1:2], ncol=2)])
      }
      
      l_M_comp[[i]] <- M_comp
      
    } # end for: i in 1:p
    
    M_comp_cmb <- do.call(rbind, l_M_comp)
    M_comp_cmb[, 1:2] <- t(apply(M_comp_cmb, 1, function(x) sort(x[1:2]) ))
    M_comp_cmb <- M_comp_cmb[!duplicated(matrix(M_comp_cmb[, 1:2], ncol=2)), ]
    M_comp_cmb <- matrix(M_comp_cmb, ncol=8)
    n_sig <- nrow(M_comp_cmb)
    
    if(n_sig == 0) {
      
      a_res[iter, 3, 5, n, 1] <- 0
      a_res[iter, 4, 5, n, 1] <- 0
      a_res[iter, 3, 5, n, 2] <- NA
      a_res[iter, 4, 5, n, 2] <- NA
      
      
    } else {
      
      # Sensitivity for Mod (with Pw)
      M_comp_cmb[, 5:6] <- matrix(rep(int_type_2[1, ], n_sig), ncol=2, byrow = TRUE)
      M_comp_cmb[, 7:8] <- matrix(rep(int_type_2[2, ], n_sig), ncol=2, byrow = TRUE)
      v_pwMod_sen <- apply(M_comp_cmb, 1, function(x) x[4] == 1 & (identical(x[1:2], x[5:6]) | identical(x[1:2], x[7:8])) )
      
      # Sensitivty for Mod (alone)
      M_comp_cmb[, 5:6] <- matrix(rep(int_type_3[1, ], n_sig), ncol=2, byrow = TRUE)
      M_comp_cmb[, 7:8] <- matrix(rep(int_type_3[2, ], n_sig), ncol=2, byrow = TRUE)
      v_Mod_sen <- apply(M_comp_cmb, 1, function(x) x[4] == 1 & (identical(x[1:2], x[5:6]) | identical(x[1:2], x[7:8])) )
      
      
      a_res[iter, 3, 5, n, 1] <- sum(v_pwMod_sen) / 2 
      a_res[iter, 4, 5, n, 1] <- sum(v_Mod_sen) / 2
      a_res[iter, 3, 5, n, 2] <- (sum(v_pwMod_sen) + sum(v_Mod_sen)) / n_sig
      a_res[iter, 4, 5, n, 2] <- (sum(v_pwMod_sen) + sum(v_Mod_sen)) / n_sig
      
    }
    
    # ----- FGL (1); moderators = m ------    
    
    M <- l_est[[iter]][[n]]$FGL_1
    m_diff <- M$part1 - M$part0
    m_diff[upper.tri(m_diff)] <- 0 # avoid getting edges twice
    ind_diff <- which(m_diff != 0, arr.ind = TRUE)
    if(is.null(dim(ind_diff))) n_diffs <- 0 else n_diffs <- nrow(ind_diff)
    
    if(n_diffs == 0) {
      
      a_res[iter, 3, 6, n, 1] <- 0
      a_res[iter, 4, 6, n, 1] <- 0
      a_res[iter, 3, 6, n, 2] <- NA
      a_res[iter, 4, 6, n, 2] <- NA
      
    } else {
      
      M_comp <- matrix(NA, nrow = n_diffs, ncol = 7)
      M_comp[, 1:2] <- ind_diff
      M_comp[, 3] <- m_diff[ind_diff]
      M_comp[, 1:2] <- t(apply(matrix(M_comp[, 1:2], ncol=2), 1, function(x) sort(x)))
      
      # Sensitivity for Mod (with Pw)
      M_comp[, 4:5] <- matrix(rep(int_type_2[1, ], n_diffs), ncol=2, byrow = TRUE)
      M_comp[, 6:7] <- matrix(rep(int_type_2[2, ], n_diffs), ncol=2, byrow = TRUE)
      v_pwMod_sen <- apply(M_comp, 1, function(x) x[3] > 0 & (identical(x[1:2], x[4:5]) | identical(x[1:2], x[6:7])) )
      
      # Sensitivty for Mod (alone)
      M_comp[, 4:5] <- matrix(rep(int_type_3[1, ], n_diffs), ncol=2, byrow = TRUE)
      M_comp[, 6:7] <- matrix(rep(int_type_3[2, ], n_diffs), ncol=2, byrow = TRUE)
      v_Mod_sen <- apply(M_comp, 1, function(x) x[3] > 0 & (identical(x[1:2], x[4:5]) | identical(x[1:2], x[6:7])) )
      
      a_res[iter, 3, 6, n, 1] <- sum(v_pwMod_sen) / 2 
      a_res[iter, 4, 6, n, 1] <- sum(v_Mod_sen) / 2
      a_res[iter, 3, 6, n, 2] <- (sum(v_pwMod_sen) + sum(v_Mod_sen)) / n_diffs
      a_res[iter, 4, 6, n, 2] <- (sum(v_pwMod_sen) + sum(v_Mod_sen)) / n_diffs
      
    }
    
    
    # ----- FGL (2); moderators = 1:p, sequentially ------
    
    # Combine Difference Estimates
    
    l_M <- l_est[[iter]][[n]]$FGL_2
    l_M_comp <- list()
    
    for(i in 1:p) {
      
      M <- l_M[[i]]
      m_diff <- M$part1 - M$part0
      m_diff[upper.tri(m_diff)] <- 0 # avoid getting edges twice
      ind_diff <- which(m_diff != 0, arr.ind = TRUE)
      n_diffs <- nrow(ind_diff)
      
      M_comp <- matrix(NA, nrow = n_diffs, ncol = 7)
      
      if(n_diffs != 0) {
        
        M_comp[, 1:2] <- ind_diff
        M_comp[, 3] <- m_diff[ind_diff]
        M_comp[, 1:2] <- t(apply(matrix(M_comp[, 1:2], ncol=2), 1, function(x) sort(x)))  
        
      } 
      
      l_M_comp[[i]] <- M_comp
      
    } # end for: i in 1:p
    
    M_comp_cmb <- do.call(rbind, l_M_comp)
    M_comp_cmb <- matrix(M_comp_cmb, ncol=7, byrow=FALSE)
    M_comp_cmb[, 1:2] <- t(apply(matrix(M_comp_cmb[, 1:2], ncol=2, byrow=FALSE), 1, sort))
    
    # order by absolute value of difference! (otherwise I overwrite correct difference estimate with incorrect difference estimate)
    M_comp_cmb <- matrix(M_comp_cmb[order(abs(M_comp_cmb[, 3]), decreasing = TRUE), ], ncol=7, byrow=FALSE)
    
    M_comp_cmb <- M_comp_cmb[!duplicated(matrix(M_comp_cmb[, 1:2], ncol=2, byrow = FALSE)), ]
    M_comp_cmb <- matrix(M_comp_cmb, ncol=7, byrow=FALSE)
    n_diffs <- nrow(M_comp_cmb) # only need n_diffs for calculation of precision
    
    
    if(n_diffs == 0) {
      
      a_res[iter, 3, 7, n, 1] <- 0
      a_res[iter, 4, 7, n, 1] <- 0
      a_res[iter, 3, 7, n, 2] <- NA
      a_res[iter, 4, 7, n, 2] <- NA
      
    } else {
      
      # Sensitivity for Mod (with Pw)
      M_comp_cmb[, 4:5] <- matrix(rep(int_type_2[1, ], n_diffs), ncol=2, byrow = TRUE)
      M_comp_cmb[, 6:7] <- matrix(rep(int_type_2[2, ], n_diffs), ncol=2, byrow = TRUE)
      v_pwMod_sen <- apply(M_comp_cmb, 1, function(x) x[3] > 0 & (identical(x[1:2], x[4:5]) | identical(x[1:2], x[6:7])) )
      
      # Sensitivty for Mod (alone)
      M_comp_cmb[, 4:5] <- matrix(rep(int_type_3[1, ], n_diffs), ncol=2, byrow = TRUE)
      M_comp_cmb[, 6:7] <- matrix(rep(int_type_3[2, ], n_diffs), ncol=2, byrow = TRUE)
      v_Mod_sen <- apply(M_comp_cmb, 1, function(x) x[3] > 0 & (identical(x[1:2], x[4:5]) | identical(x[1:2], x[6:7])) )
      
      
      a_res[iter, 3, 7, n, 1] <- sum(v_pwMod_sen) / 2 
      a_res[iter, 4, 7, n, 1] <- sum(v_Mod_sen) / 2
      a_res[iter, 3, 7, n, 2] <- (sum(v_pwMod_sen) + sum(v_Mod_sen)) / n_diffs
      a_res[iter, 4, 7, n, 2] <- (sum(v_pwMod_sen) + sum(v_Mod_sen)) / n_diffs
      
    }
    
  } # end for: n
  # print(iter)
} # end for: iter


saveRDS(a_res, file="data/ResultsArray.RDS")
a_res <- readRDS(file="data/ResultsArray.RDS")


# -----------------------------------------------------------------
# -----------------------------------------------------------------
# ----------------- Plotting --------------------------------------
# -----------------------------------------------------------------
# -----------------------------------------------------------------

# -----------------------------------------------------------------
# ----------------- Main Resuls figure ----------------------------
# -----------------------------------------------------------------



# Global Settings
cex.text <- 1.5
n_seq <- round(exp(seq(3.4, 7.5, length=12))) # n-variations
thresh_nzp <- .95 # if more than 95% are NA


# Compute Means
a_res_mean <- apply(a_res, 2:5, function(x) mean(x, na.rm = TRUE))
a_res_NAprop <- apply(a_res, 2:5, function(x) mean(is.na(x)))
a_res_NAprop[a_res_NAprop > thresh_nzp] <- NA
a_res_NAprop[a_res_NAprop <= thresh_nzp] <- 0
a_res_mean <- a_res_mean + a_res_NAprop


## Plot PDF

# Figure 5 & 6 in paper

pdf("figures/Fig3_Simulation_Results.pdf", width = 6.4, height = 8)

# Layout
lmat <- matrix(NA, nrow=5, ncol=3)
lmat[1, ] <- 1:3
lmat[2:5, 1] <- 4:7 
lmat[2:5, 2:3] <- (1:8)+7

ls <- layout(mat = lmat, 
             widths = c(.8, 1, 1), 
             heights = c(.3, 1, 1, 1, 1))

# layout.show(ls)


# Plotting

# A) Headers

par(mar = rep(0,4))
plot.new() # Empty

plot.new()
plot.window(xlim=c(-1, 1), ylim = c(-1, 1))
text(0,0, "    Sensitivity", cex = cex.text)

plot.new()
plot.window(xlim=c(-1, 1), ylim = c(-1, 1))
text(0,0, "    Precision", cex = cex.text)


# B) Left column

plot.new()
plot.window(xlim=c(-1, 1), ylim = c(-1, 1))
# text(0,0, "PW", cex = cex.text, srt=90)

plot.new()
plot.window(xlim=c(-1, 1), ylim = c(-1, 1))
# text(0,0, "PW (with Mod)", cex = cex.text, srt=90)

plot.new()
plot.window(xlim=c(-1, 1), ylim = c(-1, 1))
# text(0,0, "Mod (with PW)", cex = cex.text, srt=90)

plot.new()
plot.window(xlim=c(-1, 1), ylim = c(-1, 1))
# text(0,0, "Mod", cex = cex.text, srt=90)


# C) Results

# Select Colors
library(RColorBrewer)

cols <- brewer.pal(3, "Set1")
cols[1] <- "black" # no reason to use red, better bc formulas are in red, pot. confusing

for(per in 1:2) {
  
  par(mar=c(2, 2, 0, .25))
  
  # PW
  plot.new()
  plot.window(xlim=c(1, 12), ylim = c(0, 1))
  axis(1, n_seq, at=1:12, cex = .10, labels = FALSE)
  if(per==1) axis(2, c(0, .5, 1), las=2) else axis(2, c(0, .5, 1), las=2, labels=FALSE)
  # abline(h=c(.5, 1), col="grey")
  
  lines(a_res_mean[1, 1, , per], col = cols[1])
  lines(a_res_mean[1, 2, , per], col = cols[1], lty = 2)
  lines(a_res_mean[1, 3, , per], col = cols[1], lty = 3)
  
  # PW (with mod)
  plot.new()
  plot.window(xlim=c(1, 12), ylim = c(0, 1))
  axis(1, n_seq, at=1:12, cex = .10, labels = FALSE)
  if(per==1) axis(2, c(0, .5, 1), las=2) else axis(2, c(0, .5, 1), las=2, labels=FALSE)
  
  lines(a_res_mean[2, 1, , per], col = cols[1])
  lines(a_res_mean[2, 2, , per], col = cols[1], lty = 2)
  lines(a_res_mean[2, 3, , per], col = cols[1], lty = 3)
  
  # Legend
  if(per == 2) {
    legend(3.5, .6, c("ModNW m specified",
                      "ModNW m sequential",
                      "ModNW all m",
                      "NCT m specified",
                      "NCT m sequential",
                      "FGL m specified",
                      "FGL m sequential"),
           lty=c(1,2,3,1,2,1,2), col=c(cols[1], cols[1], cols[1], cols[2], cols[2], cols[3], cols[3]), 
           cex = .75, bty = "n")
  }
  
  # Mod (with PW)
  plot.new()
  plot.window(xlim=c(1, 12), ylim = c(0, 1))
  axis(1, n_seq, at=1:12, cex = .10, labels = FALSE)
  if(per==1) axis(2, c(0, .5, 1), las=2) else axis(2, c(0, .5, 1), las=2, labels=FALSE)
  
  lines(a_res_mean[3, 1, , per], col = cols[1])
  lines(a_res_mean[3, 2, , per], col = cols[1], lty = 2)
  lines(a_res_mean[3, 3, , per], col = cols[1], lty = 3)
  lines(a_res_mean[3, 4, , per], col = cols[2], lty = 1)
  lines(a_res_mean[3, 5, , per], col = cols[2], lty = 2)
  lines(a_res_mean[3, 6, , per], col = cols[3], lty = 1)
  lines(a_res_mean[3, 7, , per], col = cols[3], lty = 2)
  
  
  # Mod (without PW)
  plot.new()
  plot.window(xlim=c(1, 12), ylim = c(0, 1))
  axis(1, n_seq, at=1:12, cex = .10, labels = FALSE)
  if(per==1) axis(2, c(0, .5, 1), las=2) else axis(2, c(0, .5, 1), las=2, labels=FALSE)
  axis(1, n_seq, at=1:12, cex.axis = .4)
  
  lines(a_res_mean[4, 1, , per], col = cols[1])
  lines(a_res_mean[4, 2, , per], col = cols[1], lty = 2)
  lines(a_res_mean[4, 3, , per], col = cols[1], lty = 3)
  lines(a_res_mean[4, 4, , per], col = cols[2], lty = 1)
  lines(a_res_mean[4, 5, , per], col = cols[2], lty = 2)
  lines(a_res_mean[4, 6, , per], col = cols[3], lty = 1)
  lines(a_res_mean[4, 7, , per], col = cols[3], lty = 2)
  
}


dev.off()



# -----------------------------------------------------------------
# ----------------- ModNW: Sensitivity Across Par Types -----------
# -----------------------------------------------------------------

# Figure 11 in paper (Appendix)

library(RColorBrewer)

cols <- brewer.pal(4, "Set1")

# --------- Mean Sensitivity ----------

sc <- .75
lwd <- 1.5

pdf(paste0("figures/ModNW_Sensitivity_AcrossTypes.pdf"), width = sc*12, height = sc*4)

a_res_mean <- apply(a_res, 2:5, function(x) mean(x, na.rm = TRUE))

per <- 1
par(mfrow=c(1,3))

estimators <- c("(a) m specified", "(b) m sequential", "(c) all m simultaneous")

for(est in 1:3) {
  if(est ==1) ylab="Sensitivity" else ylab = ""
  
  
  plot.new()
  plot.window(ylim=c(0,1), xlim=c(1,12))
  
  axis(1, n_seq, at=1:12, cex.axis=.5)
  axis(2, c(0, .5, 1), las=2)

  lines(1:12, a_res_mean[1, est, , per], col=cols[1], lty=1, lwd=lwd)
  lines(1:12, a_res_mean[2, est, , per], col=cols[2], lty=2, lwd=lwd)
  lines(1:12, a_res_mean[3, est, , per], col=cols[3], lty=3, lwd=lwd)
  lines(1:12, a_res_mean[4, est, , per], col=cols[4], lty=4, lwd=lwd)
  
  title(ylab = ylab)
  title(xlab = "Number of observations")
  title(main = estimators[est])
  
  if(est==1) legend("topleft", c("PW", "PW(mod)", "Mod(pw)", "Mod"), 
                    lty = 1:4,
                    col = cols,
                    bty="n", 
                    lwd = rep(lwd, 4))
  
}

dev.off()

