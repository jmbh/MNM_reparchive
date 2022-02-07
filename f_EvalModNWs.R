
f_EvalModNWs <- function(est_2way, 
                         est_3way, 
                         int_type_1,
                         int_type_2,
                         int_type_3, 
                         m = m)
{
  
  
  outlist <- list("type1_S" = NULL,
                  "type1_P" = NULL,
                  "type2_S" = NULL,
                  "type2_P" = NULL,
                  "type3_S" = NULL,
                  "type3_P"= NULL,
                  "type4_S" = NULL,
                  "type4_P"= NULL)
  
  # browser()
  
  # Sensitivity: Pairwise (alone)
  M <- matrix(est_2way, ncol=2)
  n_est <- nrow(M)
  
  M_comp <- matrix(NA, nrow = n_est, ncol = 6)
  M_comp[, 1:2] <- M
  M_comp[, 3:4] <- matrix(rep(int_type_1[1, ], n_est), ncol=2, byrow = TRUE)
  M_comp[, 5:6] <- matrix(rep(int_type_1[2, ], n_est), ncol=2, byrow = TRUE)
  
  v_pw_sen <- apply(M_comp, 1, function(x) identical(x[1:2], x[3:4]) | identical(x[1:2], x[5:6]) )
  outlist$type1_S <- sum(v_pw_sen) / 2 # number of true PW=2
  
  # Sensitivity: Pairwise (in moderation)
  M_comp <- matrix(NA, nrow = n_est, ncol = 6)
  M_comp[, 1:2] <- M
  M_comp[, 3:4] <- matrix(rep(int_type_2[1, ], n_est), ncol=2, byrow = TRUE)
  M_comp[, 5:6] <- matrix(rep(int_type_2[2, ], n_est), ncol=2, byrow = TRUE)
  
  v_pwMod_sen <- apply(M_comp, 1, function(x) identical(x[1:2], x[3:4]) | identical(x[1:2], x[5:6]) )
  outlist$type2_S <- sum(v_pwMod_sen) / 2 # number of true PW=2
  
  # Precision: Pairwise (both alone and in moderation; the same)
  outlist$type1_P <- (sum(v_pw_sen) + sum(v_pwMod_sen)) / n_est
  outlist$type2_P <- (sum(v_pw_sen) + sum(v_pwMod_sen)) / n_est
  
  # Sensitivity: Moderator (alone)
  M <- matrix(est_3way, ncol = 3)
  M <- matrix(apply(M, 1, function(x) sort(x)), ncol=3, byrow = T) # sort
  
  n_est <- nrow(M)
  
  M_comp <- matrix(NA, nrow = n_est, ncol = 9)
  M_comp[, 1:3] <- M
  fill_A <- matrix(rep(c(int_type_2[1, ], m), n_est), ncol=3, byrow = TRUE)
  fill_A <- matrix(apply(fill_A, 1, function(x) sort(x)), ncol=3, byrow = TRUE)
  M_comp[, 4:6] <- fill_A
  fill_B <- matrix(rep(c(int_type_2[2, ], m), n_est), ncol=3, byrow = TRUE)
  fill_B <- matrix(apply(fill_B, 1, function(x) sort(x)), ncol=3, byrow = TRUE)
  M_comp[, 7:9] <- fill_B
  
  v_Mod_sen <- apply(M_comp, 1, function(x) identical(x[1:3], x[4:6]) | identical(x[1:3], x[7:9]) )
  outlist$type3_S <- sum(v_Mod_sen) / 2 # number of true PW=2
  
  
  # Sensitivity: Moderator (in pairwise)
  M_comp <- matrix(NA, nrow = n_est, ncol = 9)
  M_comp[, 1:3] <- M
  fill_A <- matrix(rep(c(int_type_3[1, ], m), n_est), ncol=3, byrow = TRUE)
  fill_A <- matrix(apply(fill_A, 1, function(x) sort(x)), ncol=3, byrow = TRUE)
  M_comp[, 4:6] <- fill_A
  fill_B <- matrix(rep(c(int_type_3[2, ], m), n_est), ncol=3, byrow = TRUE)
  fill_B <- matrix(apply(fill_B, 1, function(x) sort(x)), ncol=3, byrow = TRUE)
  M_comp[, 7:9] <- fill_B
  
  
  v_ModPW_sen <- apply(M_comp, 1, function(x) identical(x[1:3], x[4:6]) | identical(x[1:3], x[7:9]) )
  outlist$type4_S <- sum(v_ModPW_sen) / 2 # number of true PW=2
  
  # Precision: Moderator (both alone and in moderation; the same)
  outlist$type3_P  <- (sum(v_Mod_sen) + sum(v_ModPW_sen)) / n_est
  outlist$type4_P  <- (sum(v_Mod_sen) + sum(v_ModPW_sen)) / n_est
  
  
  return(outlist)
  
  
} # end function: f_EvalModNWs
