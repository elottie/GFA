gfa_singletrait_check <- function(fit, check_thresh = 0.9, params){

  D <- fit$F_pm
  nfactors <- ncol(D)
  L <- fit$L_pm

  if(length(fit$flash_fit$fix.dim) == 0){
    fixed_ix <- rep(FALSE, nfactors)
  }else{
    fixed_ix <- fit$flash_fit$fix.dim %>% sapply(., function(x){
      if(is.null(x)) return(FALSE)
      if(x == dim) return(TRUE)
      return(FALSE)})
  }

  Dn <- norm_cols(D)$A
  col_max <- apply(abs(Dn), 2, max)

  single_trait_index <- which(col_max > check_thresh & !fixed_ix)

  for(i in single_trait_index){
    cat("Checking factor ", i, "\n")
    myfactor <- D[,i, drop = FALSE]
    myloadings <- L[, i, drop = FALSE]
    altfactor <- myfactor
    k <- which.max(abs(myfactor))
    altfactor[-k] <- 0

    new_order <- seq(ncol(D))
    if(i < nfactors){
      new_order[i:(nfactors-1)] <- (i:(nfactors-1)) + 1
      new_order[nfactors] <- i
    }
    fitn <- flash_factors_remove(fit, i) %>%
            flash_factors_init(init = list(myloadings, altfactor),
                               ebnm_fn = list(params$ebnm_fn_L, params$ebnm_fn_F)) %>%
            #flash_factors_reorder(new_order) %>%
            flash_factors_fix(., kset = nfactors, which_dim = "factors") %>%
            flash_backfit()

    if(i < nfactors){
      fitn <- flash_factors_reorder(fitn, new_order)
    }
    if(fitn$elbo > fit$elbo){
      message(paste0("Replacing factor ", i , " with single trait factor"))
      fit <- fitn
    }
  }

  return(fit)
}
