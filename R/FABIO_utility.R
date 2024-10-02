# function to prepare GReX for a single gene
prepGReX_gene = function(i){
  tb_weight = read.table(paste0(weight_dir,'/',weight_list[i]))
  gene = unique(tb_weight[,1])
  pos = strsplit(tb_weight[,2],':')
  pos = sapply(pos,'[[',2)
  pos = as.numeric(pos)
  end1 = min(pos)
  end2 = max(pos)
  
  # load genotype of selected chr and region
  plink_cmd = paste0("plink2 --bfile ",geno_dir,"/chr",chr," --chr ",chr," --from-bp ",end1," --to-bp ",end2," --make-bed --out ./temp_files/",gene)
  system(plink_cmd)
  
  geno_file = read.plink(paste0("./temp_files/",gene,".bed"))
  snp_id = paste0(geno_file$map$chromosome,":",geno_file$map$position)
  genotype_df = t(scale(as(geno_file$genotypes, "numeric")))
  genotype_df[is.na(genotype_df)] = 0
  rownames(genotype_df) = snp_id
  genotype_df = genotype_df[tb_weight[,2],]
  
  grex = as.vector(t(tb_weight[,5])%*%genotype_df)
  return(c(gene,grex))
}

# function to calculate Pgamma
CalcPgamma = function(geo_mean, ng_test) {
  p_gamma = 0:(ng_test - 1)
  p_gamma = 0.7 * dgeom(p_gamma, 1.0 / geo_mean, log = FALSE) + 0.3 / ng_test
  s = sum(p_gamma)
  p_gamma = p_gamma / s
  return(p_gamma)
}

# function to initiate MCMC
InitialMCMC = function(pos_loglr, ng_test, g_max, g_min, pve_null, h_min, h_max, logp_min, logp_max) {
  q_genome = qchisq(0.05 / ng_test, 1, lower.tail = F)
  cHyp = list(n_gamma = 0)
  for (i in 1:nrow(pos_loglr)) {
    if (2 * pos_loglr[i,2] > q_genome) {
      cHyp$n_gamma = cHyp$n_gamma + 1
    }
  }
  if (cHyp$n_gamma < 10) {cHyp$n_gamma = 10}
  if (cHyp$n_gamma > g_max) {cHyp$n_gamma = g_max}
  if (cHyp$n_gamma < g_min) {cHyp$n_gamma = g_min}

  rank = 1:cHyp$n_gamma

  cHyp$logp = log(cHyp$n_gamma / ng_test)
  cHyp$h = pve_null

  if (cHyp$logp == 0) {cHyp$logp = -0.000001}
  if (cHyp$h == 0) {cHyp$h = 0.1}

  if (cHyp$h < h_min) {cHyp$h = h_min}
  if (cHyp$h > h_max) {cHyp$h = h_max}

  if (cHyp$logp < logp_min) {cHyp$logp = logp_min}
  if (cHyp$logp > logp_max) {cHyp$logp = logp_max}

  cHyp$pve = 0
  cHyp$pge = 0
  return(list(cHyp=cHyp, rank=rank))
}

# function to propose the new gamma
ProposeGamma = function(rank_old, p_gamma, cHyp_old, cHyp_new, rep, g_max, g_min, ng_test) {
  # sample rv from discrete
  sample_discrete = function(n) {
    sample(x = 1:ng_test, n, replace = T, prob = p_gamma)
  }

  unif = 0
  logp = 0

  mapRank2in = data.frame(rank_old=integer(), check=integer())
  rank_new = c()

  if (cHyp_old$n_gamma != 0) {
    rank_new = rep(NA, length(rank_old))
    mapRank2in = rbind(mapRank2in, cbind(rank_old, NA))
    for (i in 1:length(rank_old)) {
      r = rank_old[i]
      rank_new[i] = r
      mapRank2in[i,2] = 1
    }
  }
  cHyp_new$n_gamma = cHyp_old$n_gamma

  for (i in 1:rep) {
    unif = runif(1)

    if (unif < 0.40 && cHyp_new$n_gamma < g_max) {
      flag_gamma = 1
    } else if (unif >= 0.40 && unif < 0.80 && cHyp_new$n_gamma > (g_min+2)) {
      flag_gamma = 2
    } else if (unif >= 0.80 && cHyp_new$n_gamma > 0 && cHyp_new$n_gamma < ng_test) {
      flag_gamma = 3
    } else {
      flag_gamma = 4
    }

    if (flag_gamma == 1) {
      # add a SNP
      r_add = sample_discrete(1)
      while (r_add %in% mapRank2in[,1] == T) {r_add = sample_discrete(1)}

      prob_total = 1
      if (cHyp_new$n_gamma != 0){
        for (i in 1:cHyp_new$n_gamma) {
          r = rank_new[i]
          prob_total = prob_total - p_gamma[r]
        }
      }

      mapRank2in = rbind(mapRank2in, c(r_add,1))
      rank_new = append(rank_new, r_add)
      cHyp_new$n_gamma = cHyp_new$n_gamma + 1
      logp = logp - log(p_gamma[r_add] / prob_total) - log(cHyp_new$n_gamma)
    } else if (flag_gamma == 2) {
      # Delete a SNP
      col_id = 1 + floor(runif(1, 0, cHyp_new$n_gamma))
      r_remove = rank_new[col_id]

      prob_total = 1
      for (i in 1:cHyp_new$n_gamma) {
        r = rank_new[i]
        prob_total = prob_total - p_gamma[r]
      }
      prob_total = prob_total + p_gamma[r_remove]

      if (cHyp_new$n_gamma == 1) {
        mapRank2in = data.frame(rank_old=integer(), check=integer())
      } else {
        mapRank2in = mapRank2in[!mapRank2in[,1]==r_remove, ]
      }
      rank_new = setdiff(rank_new, rank_new[col_id])
      if(length(rank_new)==0) {rank_new = c()}

      logp = logp + log(p_gamma[r_remove] / prob_total) + log(cHyp_new$n_gamma)
      cHyp_new$n_gamma = cHyp_new$n_gamma - 1
    } else if (flag_gamma == 3) {
      # Switch a SNP
      col_id = 1 + floor(runif(1, 0, cHyp_new$n_gamma))
      r_remove = rank_new[col_id]

      r_add = sample_discrete(1)
      while (r_add %in% mapRank2in[,1] == T) {r_add = sample_discrete(1)}

      prob_total = 1
      for (i in 1:cHyp_new$n_gamma) {
        r = rank_new[i]
        prob_total = prob_total - p_gamma[r]
      }

      logp = logp + log(p_gamma[r_remove] / (prob_total + p_gamma[r_remove] - p_gamma[r_add]))
      logp = logp - log(p_gamma[r_add] / prob_total)

      if (cHyp_new$n_gamma == 1) {
        mapRank2in = data.frame(rank_old=integer(), check=integer())
      } else {
        mapRank2in = mapRank2in[!mapRank2in[,1]==r_remove, ]
      }
      mapRank2in = rbind(mapRank2in, c(r_add,1))
      rank_new = setdiff(rank_new, rank_new[col_id])
      if(length(rank_new)==0) {rank_new = c()}
      rank_new = append(rank_new, r_add)
    } else {
      logp = logp + 0
    }
  }

  rank_new = sort(rank_new)
  return(list(logp=logp, rank_new=rank_new, cHyp_new=cHyp_new))
}

# function to set up Xgamma matrix
SetXgamma = function(X, X_old, XtX_old, Xty_old, y, rank_old, rank_new, X_new, XtX_new, Xty_new) {
  # rank_old and rank_new are sorted already inside ProposeGamma
  # calculate vectors rank_remove and rank_add.
  rank_remove = setdiff(rank_old, rank_new)
  if (identical(rank_remove, integer(0))) {rank_remove = 0}
  rank_add = setdiff(rank_new, rank_old)
  if (identical(rank_add, integer(0))) {rank_add = 0}
  rank_union = union(rank_new, rank_old)

  # Map rank_remove and rank_add.
  mapRank2in_remove = cbind(rank_remove, 1)
  if (length(rank_remove) == 1) {
    if (rank_remove == 0) {rank_remove = c()}
  }
  mapRank2in_add = cbind(rank_add, 1)
  if (length(rank_add) == 1) {
    if (rank_add == 0) {rank_add = c()}
  }

  # Obtain the subset of matrix/vector
  Xold_sub = X_old[,1:length(rank_old)]
  XtXold_sub = XtX_old[1:length(rank_old), 1:length(rank_old)]
  Xtyold_sub = Xty_old[1:length(rank_old)]

  Xnew_sub = X_new[,1:length(rank_new)]
  XtXnew_sub = XtX_new[1:length(rank_new), 1:length(rank_new)]
  Xtynew_sub = Xty_new[1:length(rank_new)]

  # Get X_new and calculate XtX_new
  if (length(rank_remove) == 0 & length(rank_add) == 0) {
    X_new[,1:length(rank_new)] = Xold_sub
    XtX_new[1:length(rank_new), 1:length(rank_new)] = XtXold_sub
    Xty_new[1:length(rank_new)] = Xtyold_sub
  } else {
    if (length(rank_add) == 0) {
      i_old = 1
      i_new = 1
      for (i in 1:length(rank_union)) {
        if (rank_old[i_old] %in% mapRank2in_remove[,1]) {
          i_old = i_old + 1
          next
        }

        X_new[, i_new] = X_old[, i_old]

        Xty_new[i_new] = Xty_old[i_old]

        j_old = i_old
        j_new = i_new
        for (j in i:length(rank_union)) {
          if (rank_old[j_old] %in% mapRank2in_remove[,1]) {
            j_old = j_old + 1
            next
          }

          d = XtX_old[i_old, j_old]
          XtX_new[i_new, j_new] = d
          if (i_new != j_new) {XtX_new[j_new, i_new] = d}

          j_old = j_old + 1
          j_new = j_new + 1
        }
        i_old = i_old + 1
        i_new = i_new + 1
      }
    } else {
      X_add = matrix(NA, nrow(X_old), length(rank_add))
      XtX_aa = matrix(NA, ncol(X_add), ncol(X_add))
      XtX_ao = matrix(NA, ncol(X_add), ncol(X_old))
      Xty_add = rep(NA, ncol(X_add))
      # Get X_add.
      X_add[,1:length(rank_add)] = X[,mapRank2pos[rank_add,2]]

      XtX_aa = t(X_add) %*% X_add
      XtX_ao = t(X_add) %*% X_old
      Xty_add = t(X_add) %*% y

      # Save to X_new, XtX_new and Xty_new.
      i_old = 1
      i_new = 1
      i_add = 1
      for (i in 1: length(rank_union)) {
        if (rank_old[i_old] %in% mapRank2in_remove[,1]) {
          i_old = i_old + 1
          next
        }
        if (rank_new[i_new] %in% mapRank2in_add[,1]) {
          i_flag = 1
        } else {
          i_flag = 0
        }

        if (i_flag == 1) {
          Xcopy_col = X_add[, i_add]
          X_new[, i_new] = Xcopy_col
        } else {
          Xcopy_col = X_old[, i_old]
          X_new[, i_new] = Xcopy_col
        }

        if (i_flag == 1) {
          d = Xty_add[i_add]
        } else {
          d = Xty_old[i_old]
        }
        Xty_new[i_new] = d

        j_old = i_old
        j_new = i_new
        j_add = i_add
        for (j in i:length(rank_union)) {
          if (rank_old[j_old] %in% mapRank2in_remove[,1]) {
            j_old = j_old + 1
            next
          }
          if (rank_new[j_new] %in% mapRank2in_add[,1]) {
            j_flag = 1
          } else {
            j_flag = 0
          }

          if (i_flag == 1 & j_flag == 1) {
            d = XtX_aa[i_add, j_add]
          } else if (i_flag == 1) {
            d = XtX_ao[i_add, j_old]
          } else if (j_flag == 1) {
            d = XtX_ao[j_add, i_old]
          } else {
            d = XtX_old[i_old, j_old]
          }

          XtX_new[i_new, j_new] = d
          if (i_new != j_new) {
            XtX_new[j_new, i_new] = d
          }

          j_new = j_new + 1
          if (j_flag == 1) {
            j_add = j_add + 1
          } else {
            j_old = j_old + 1
          }
        }
        i_new = i_new + 1
        if (i_flag == 1) {
          i_add = i_add + 1
        } else {
          i_old = i_old + 1
        }
      }
    }
  }
  return(list(X_new=X_new, XtX_new=XtX_new, Xty_new=Xty_new))
}
