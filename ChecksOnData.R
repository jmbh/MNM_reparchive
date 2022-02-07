# jonashaslbeck@gmail.com; February 2018 


# -----------------------------------------------------------------
# ----------------- Load Results ----------------------------------
# -----------------------------------------------------------------


# !!! Specity Directories of Simulation Output !!!

outdir_data <- ".../output_data/" # Specify directory of data output of simulation
outdir_est <- ".../output_est/" # Specify directory of estimation output of simulation


# ----- Load Data files -----

dfiles_data <- list.files(outdir_data)
n_files_data <- length(dfiles_data)

# get iters_data
iters_data <- substr(dfiles_data, 18, 20)
iters_data <- sub(".", "", iters_data, fixed=TRUE)
iters_data <- sub("R", "", iters_data, fixed=TRUE)
iters_data <- as.numeric(iters_data)


# ----- Load Estimation files -----

dfiles_est <- list.files(outdir_est)
n_files_est <- length(dfiles_est)

# get iters_est
iters_est <- substr(dfiles_est, 17, 19)
iters_est <- sub(".", "", iters_est, fixed=TRUE)
iters_est <- sub("R", "", iters_est, fixed=TRUE)
iters_est <- as.numeric(iters_est)


# ----- Select iters for which we have estimates -----

# Which ones are estimated?
iters_with_est <- iters_data[iters_data %in% iters_est]
dfiles_data_2load <- dfiles_data[iters_data %in% iters_est]
n_files <- length(dfiles_data_2load)

l_data <- list()
for(i in 1:n_files) l_data[[i]] <- readRDS(paste0(outdir_data, dfiles_data_2load[i]))


# -----------------------------------------------------------------
# ----------------- Get Proportion of finite for each iteration ---
# -----------------------------------------------------------------

# Get proportion of finite samples
mDiv <- matrix(NA, nrow = n_files, ncol = 2)
mDiv[, 1] <- iters_with_est
for(i in 1:n_files) mDiv[i, 2] <- l_data[[i]]$Datalist$pop_finite

ordD <- order(mDiv[, 2])
mDiv_ord <- mDiv[ordD, ]

# Select 100 iterations with largest proportion of finite samples
first100 <- mDiv_ord[-(1:33) ,] # take 100
v_first100 <- first100[, 1]
saveRDS(v_first100, file = "data/first100.RDS") # Will be usd in modSim_eval.R


# -----------------------------------------------------------------
# ----------------- Plot all 100 used Graphs  ---------------------
# -----------------------------------------------------------------

# This is not shown in the paper; but might be interesting to inespect

# look at graphs of extremes
library(qgraph)

ind_display <- cbind(1:n_files, iters_with_est, mDiv[, 2])
ind_display <- ind_display[iters_with_est %in% v_first100, ]
ind_display <- ind_display[order(ind_display[,3]), ]

pdf("figures/ModSim_PlotAllGraphs.pdf", width = 8 *2, height = 7*2)

par(mfrow=c(10,10))

for(i in 1:100) {
  
  iter <- ind_display[i, 2]
  file_i <- ind_display[i, 1]
  
  
  node.color <- rep("white", 13)
  node.color[l_data[[file_i]]$m] <- "grey"
  
  edge.color <- matrix("green", 13, 13)
  edge.color[l_data[[file_i]]$Graph == 2] <- "orange"
  edge.color[l_data[[file_i]]$Graph == 3] <- "blue"
  
  qgraph(l_data[[file_i]]$Graph, 
         weighted = FALSE,
         edge.labels = TRUE, 
         title = paste0("Iteration = ", iter, ", Convergence = ", ind_display[i, 3]),
         title.cex = .7,
         color = node.color, 
         edge.color = edge.color,
         edge.label.cex = 4,
         layout = "circle")
  
}

dev.off()


# -----------------------------------------------------------------
# ----------------- Check on Bias ---------------------------------
# -----------------------------------------------------------------

# Reload: now only the 100 used
dfiles_data_2load <- dfiles_data[iters_data %in% v_first100]
n_files <- length(dfiles_data_2load)

l_data <- list()
for(i in 1:n_files) l_data[[i]] <- readRDS(paste0(outdir_data, dfiles_data_2load[i]))

n_seq <- round(exp(seq(3.4, 7.5, length=12))) # n-variations

for(n in 12) {
  
  # iterations, interaction types, estimates
  a_bias <- array(NA, dim = c(n_files, 4, 2))
  
  for(i in 1:n_files) {
    
    data <- l_data[[i]]$Datalist$data_finite
    data <- data[1:n_seq[n], ]
    G <- l_data[[i]]$Graph
    G[upper.tri(G)] <- 0
    m <- l_data[[i]]$m
    
    data <- as.data.frame(data)
    
    counter_t1 <- 1
    counter_t2 <- 1
    counter_t3 <- 1
    
    for(v in 1:13) {
      
      fitobj <- lm(as.formula(paste0(colnames(data)[v], "~ (.)^2")), data = data)  
      C <- fitobj$coefficients
      
      ind_inters <- which(G[v, ] != 0)
      if(length(ind_inters)>0) {
        for(a in 1:length(ind_inters)) {
          
          # type 1
          if(G[v, ][ind_inters[a]] == 1) {
            a_bias[i, 1, counter_t1] <- C[names(C) == paste0("V", ind_inters[a])] # the pairwise interaction
            counter_t1 <- counter_t1 + 1
          }
          
          # type 2
          if(G[v, ][ind_inters[a]] == 2) {
            a_bias[i, 2, counter_t2] <- C[names(C) == paste0("V", ind_inters[a])] # the pairwise interaction
            
            vars <- sort(c(ind_inters[a], m))
            a_bias[i, 3, counter_t2] <- C[names(C) ==  paste0("V", vars[1], ":V", vars[2])] # the moderation effect
            counter_t2 <- counter_t2 + 1
          }
          
          # type 3
          if(G[v, ][ind_inters[a]] == 3) {
            vars <- sort(c(ind_inters[a], m))
            a_bias[i, 4, counter_t3] <- C[names(C) ==  paste0("V", vars[1], ":V", vars[2])] # the moderation effect
            counter_t3 <- counter_t3 + 1
          }
          
          
        } # end for: interactions in row v
      } # end if: interaction in row?
      
    } # end for: rows
    
  } # end for: iter
  
  print(n)
  
  # Make Summary
  mBias <- matrix(NA, nrow = n_files, ncol = 6)
  mBias[, 1] <- v_first100
  for(i in 1:n_files) mBias[i, 2] <- l_data[[i]]$Datalist$pop_finite
  mBias[, 3:6] <- apply(a_bias, 1:2, mean)
  
  mean(mBias[, 3])
  mean(mBias[, 4])
  mean(mBias[, 5])
  mean(mBias[, 6])
  
} # end for: n


# --------- Final Figure for Paper with ideal Scaling ---------

# This is Figure 10 in the paper

# Plotting
pdf(paste0("figures/Data_SC_ModNW_n", n_seq[n], "_paperVersion.pdf"), width = 8, height = 6)

par(mfrow=c(2,2))
mains <- c("Pairwise", "Pairwise (with Moderation)", "Moderation (with Pairwise)", "Moderation alone")
for(i in 1:4) {
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0, .4))
  axis(1, seq(0, 1, length = 6))
  axis(2, c(0, .2, .4), las=2)
  
  abline(h=.2, col="black", lty=1)
  points(mBias[, 2], mBias[, i+2])
  abline(lm(mBias[, i+2]~ mBias[, 2]), lty=2, col="red")
  
  title(main = mains[i], xlab = "Proportion finite samples", ylab="Mean estimate")
}

dev.off()




