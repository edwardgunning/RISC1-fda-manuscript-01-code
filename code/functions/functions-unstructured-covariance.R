# This script comtains code from the denseFLMM package 
# (https://cran.r-project.org/web/packages/denseFLMM/index.html)
# it is re-purposed to estimate random effects covariances for mutlivariate
# functional data represented by basis expansions.

# For now, it only works for bivariate functional data with a single level
# of random intercepts, as described in the Appendix of paper.

cov_unstruct_coef <- function(Ycoef,
                                Zlist = NA, 
                                G = NA,
                                Lvec = NA,
                                groups = matrix(1, nrow(Ycoef), 1),
                                Zvars) {
  # This part is taken from the denseFLMM function.
  # -------------------------------------------------------------------------
  ###########################################################################
  # checks for consistency of input
  ###########################################################################
  if(all(is.na(Zlist)) & all(is.na(groups))){
    stop("either Zlist or groups and Zvars must be specified")
  }
  if(all(is.na(Zlist)) & all(is.na(Zvars))){
    stop("either Zlist or groups and Zvars must be specified")
  }
  if(all(is.na(Zlist))){
    if(nrow(Ycoef) != nrow(groups)){
      stop("The number of rows in Ycoef needs to agree with the number
           of rows in the grouping matrix")
    }
    if(ncol(groups) != length(Zvars)){
      stop("the number of grouping factors has to correspond to
           the number of groups of random variables")
    }
    if(!prod(sapply(seq(len = length(Zvars)),
                    function(r) nrow(Zvars[[r]])) == nrow(Ycoef))){
      stop("the number of rows in Ycoef needs to agree with the number
           of rows in the matrices of random variables")
    }
  } else{
    # each matrix in Zlist should be of dim nrow(Y) x L^U_g
    check_dims <- function(Zlist, Ycoef, Lvec){
      zdim <- list()
      for(z in seq_along(Zlist)){
        zdim[[z]] <- do.call(rbind, (lapply(Zlist[[z]], dim)))
      }
      zdim_un <- do.call(rbind, zdim)
      if(isTRUE(all.equal(zdim_un[, 1], rep(zdim_un[1, 1], nrow(zdim_un))))){
        x <- isTRUE(all.equal(zdim_un[1, 1], nrow(Ycoef)))
        if(!x)
          stop("the number of rows of each matrix in Zlist need to
                 correspond to the number of rows of Ycoef")
      }else{
        stop("the number of rows of each matrix in Zlist need to
               correspond to the number of rows of Ycoef")
      }
      y <- (all(Lvec == zdim_un[, 2]))
      if(!y)
        stop("the number of columns of each matrix in Zlist need to
               correspond to the respective number in Lvec")
    }
    check_dims(Zlist = Zlist, Y = Ycoef, Lvec = Lvec)
  }
  # if(ncol(Y) != length(gridpoints)){
  #   stop("the number of columns in Y needs to agree with the length of
  #        the gridpoints vector")
  # }
  # if(is.na(L) & (!prod(!is.na(NPC)))){
  #   warning("as both L and part of NPC are missing, will default
  #           to L = 0.9, NPC = NA")
  #   L <- 0.9
  #   NPC <- NA
  # }
  # if(!is.na(L) & (prod(!is.na(NPC)))){
  #   warning("NPC will override choice of L")
  # }
  # if(!is.na(L)){
  #   if(L > 1 | L < 0){
  #     stop("the level of explained variance needs to be between 0 and 1")
  #   }
  # }
  # if(smooth)
  #   stopifnot(smoothalg %in% c("gamm", "gamGCV", "gamREML", "bamGCV",
  #                              "bamREML", "bamfREML"))
  
  
  message("set up")
  D <- ncol(Ycoef)   # number of basis coefficients per curve
  n <- nrow(Ycoef)   # overall number of curves
  
  if(all(is.na(Zlist))){  # no group-specific functional random effects
    G <- length(Zvars)    # number of grouping factors
    rhovec <- 1
    Lvec <- n
    if(G > 0){
      # number of random effects for grouping factor 1,..., G
      rhovec <- c(sapply(1:G, function(g){ncol(Zvars[[g]])}), 1)
      
      # number of levels of grouping factors 1,..., G
      Lvec <- c(sapply(1:G, function(g){nlevels(as.factor(groups[, g]))}), n)
    }
  }else{
    H <- length(Zlist)
    rhovec <- unlist(lapply(Zlist, length))
  }
  sq2 <- sum(rhovec^2)
  
  # construct H, Zlist when no group-specific
  # functional random effects are present
  ###########################################
  if(all(is.na(Zlist))){  # define H as G + 1 (for smooth error)
    H <- G + 1
  }
  if(all(is.na(Zlist))){
    message("construct Zlist")
    foroneq <- function(g, q){
      gp1 <- as.factor(groups[, g])
      sparse.model.matrix(~ gp1 * Zvars[[g]][, q] - gp1 - Zvars[[g]][, q] - 1)
    }
    Zlist <- lapply(seq(len = G), function(g){lapply(1:rhovec[g],
                                                     function(q) foroneq(g, q))})
    Zlist[[H]] <- list(Diagonal(n))
  }
  
  # assume a centered Y for now, no estimation of mean function
  Y.tilde <- Ycoef
  
  ###########################################################################
  # estimation
  ###########################################################################
  
  # estimate covariance functions using least squares
  ###################################################
  message("estimate covariance(s)")
  gcyc <- rep(1:(H), rhovec^2)
  qcyc <- unlist(sapply(rhovec, FUN = function(q){rep(1:q, each = q)}))
  pcyc <- unlist(sapply(rhovec, FUN = function(p){rep(1:p, p)}))
  XtXentry <- function(ro, co){
    A1 <-
      if(H == G + 1){ # if no group-specific smooth errors
        # avoid unnecessary multiplications with pseudo-Design
        # Zlist[[H]] (= Identity) for residuals:
        if(gcyc[co] == H){
          Zlist[[gcyc[ro]]][[pcyc[ro]]]
        }else{
          if(gcyc[ro] == H){
            t(Zlist[[gcyc[co]]][[pcyc[co]]])
          }else{
            crossprod(Zlist[[gcyc[co]]][[pcyc[co]]],
                      Zlist[[gcyc[ro]]][[pcyc[ro]]])
          }
        }
      }else{
        crossprod(Zlist[[gcyc[co]]][[pcyc[co]]],
                  Zlist[[gcyc[ro]]][[pcyc[ro]]])
      }
    A2 <-
      if(H == G + 1){ # if no group-specific smooth errors
        # avoid unnecessary multiplications with pseudo-Design
        # Zlist[[H]] (= Identity) for residuals:
        if(gcyc[co] == H){
          Zlist[[gcyc[ro]]][[qcyc[ro]]]
        }else{
          if(gcyc[ro] == H){
            t(Zlist[[gcyc[co]]][[qcyc[co]]])
          }else{
            crossprod(Zlist[[gcyc[co]]][[qcyc[co]]],
                      Zlist[[gcyc[ro]]][[qcyc[ro]]])
          }
        }
      }else{
        crossprod(Zlist[[gcyc[co]]][[qcyc[co]]],
                  Zlist[[gcyc[ro]]][[qcyc[ro]]])
      }
    
    # get trace tr(A1'A2) without computing off-diagonal elements in A1'A2
    traceA1tA2 <- function(A1, A2){
      ret <- if(all(dim(A1) == dim(A2))){
        sum(rowSums(A1 * A2))
      }else{
        # use tr(A1'A2) = tr(A2'A1) to shorten loop
        if(ncol(A1) < nrow(A1)){
          sum(sapply(1:ncol(A1),
                     function(i) as.numeric(crossprod(A1[, i, drop = F],
                                                      A2[, i, drop = F])), simplify = T))
        }else{
          sum(sapply(1:nrow(A1),
                     function(i) as.numeric(tcrossprod(A1[i, , drop = F],
                                                       A2[i,, drop = F])), simplify = T))
        }
      }
      return(ret)
    }
    return(traceA1tA2(A1, A2))
  }
  
  # define useful function similar to outer
  # does not need function to work for vectors
  # or to return scalars
  ############################################
  matrixouter <- function(rows, cols, FUN){
    FUN <- match.fun(FUN)
    # use rbind, cbind to preserve sparsity patterns
    do.call(rbind, lapply(rows, function(g)
      do.call(cbind, lapply(cols, function(h) FUN(g, h)))))
  }
  
  XtX <- matrixouter(seq(len = sq2), seq(len = sq2), FUN = "XtXentry")
  
  Xtcentry <- function(ro){
    if(gcyc[ro] == H & H == (G + 1)){
      return(as.vector(crossprod(Y.tilde)))
    }else{
      A1 <- crossprod(Zlist[[gcyc[ro]]][[pcyc[ro]]], Y.tilde)
      A2 <- crossprod(Zlist[[gcyc[ro]]][[qcyc[ro]]], Y.tilde)
      return(as.vector(crossprod(A1, A2)))
    }
  }
  
  Xtc <- do.call(rbind, lapply(seq(len = sq2), function(h){Xtcentry(h)}))
  Ktilde <- solve(XtX, Xtc)
  # -------------------------------------------------------------------------
  # -------------------------------------------------------------------------
  onecov <- function(g){
    covcomp <- function(s, r){
      matrix(Ktilde[sum(rhovec[seq(len = g - 1)]^2) + (s - 1) * rhovec[g] + r, ],
             D, byrow = TRUE)
    }
    matrixouter(1:rhovec[g], 1:rhovec[g], FUN = "covcomp")
  }
  lapply(1:(H), "onecov")
}



# -------------------------------------------------------------------------
# For univariate case:
# -------------------------------------------------------------------------
cov_unstruc_mlfpca_uni_fd <- function(fd_obj,
                                      id_vec) {
  # going to work with the `fda` package
  require(fda) # CRAN v5.5.1
  require(Matrix) # CRAN v1.4-0
  
  # checks on the FDA object:
  stopifnot(is.fd(fd_obj))
  # check it's a univariate fda object
  stopifnot(length(dim(fd_obj$coefs)) == 2)
  
  # Grab some stuff from fd object:
  basis_obj <- fd_obj$basis
  coef_mat <- fd_obj$coef
  N <- ncol(coef_mat)
  K <- nrow(coef_mat)
  
  if (N < 2) stop("PCA not possible without replications.")
  
  # Check the id_vec:
  stopifnot(is.integer(id_vec))
  if(length(id_vec) != N) stop("Length of id_vec must match number of observations in the fd_obj")
  
  N_subject <- length(unique(id_vec))
  if(!(N_subject < N)) stop("Must have replicate observations from at least some subjects")
  
  # Create design matrix:
  id_factor <- factor(id_vec)
  Zmat <- sparse.model.matrix(~ - 1 + id_factor)
  stopifnot(dim(Zmat) == c(N, N_subject))
  
  # And put into list with indetity matrix for our covariance estimation function (requires this because based of denseFLMM)
  Z_list <- list(list(), list())
  Z_list[[1]][[1]] <- Zmat
  Z_list[[2]][[1]] <- Diagonal(n = N)
  
  # Do a least squares fit on the coefficients (Cederbaum 2017 Section 3/ Appendix A of our paper)
  K_coef_list <- cov_unstruct_coef(Ycoef = t(coef_mat),
                                   Zlist = Z_list, 
                                   G = 1, 
                                   Lvec = c(N_subject, N),
                                   groups = NA,
                                   Zvars = NA)
  
  stopifnot(length(K_coef_list) == 2)
  
  K_U_coef <- K_coef_list[[1]] # get basis coefs for the between-subject covariance function
  stopifnot(dim(K_U_coef) == rep(K, 2)) 
  K_U_bifd <- bifd(coef = K_U_coef, # and use coefficients to create bivariate fd object
                   sbasisobj = basis_obj,
                   tbasisobj = basis_obj) # i.e., Q in our paper
  
  K_E_coef <- K_coef_list[[2]] # get basis coefs for the within-subject covariance function
  stopifnot(dim(K_E_coef) == rep(K, 2)) 
  K_E_bifd <- bifd(coef = K_E_coef, 
                   sbasisobj = basis_obj,
                   tbasisobj = basis_obj) # i.e., S in our paper
  
  # Return bivariate fd object
  list(K_U_bifd = K_U_bifd, K_E_bifd = K_E_bifd)
}



# -------------------------------------------------------------------------
# For bivariate case, accepting bivariate functional data object:
# -------------------------------------------------------------------------
cov_unstruc_mlfpca_bi_fd <- function(fd_obj_list,
                                      id_vec) {
  # going to work with the `fda` package
  require(fda) # CRAN v5.5.1
  require(Matrix) # CRAN v1.4-0
  
  stopifnot(is.list(fd_obj_list)) # must be a list
  stopifnot(length(fd_obj_list) == 2) # check than our bivariate `fda` object is made up of two univariate `fda` objects
  # checks that each element of the list is an `fda` object:
  stopifnot(sapply(fd_obj_list, is.fd))
  # check each element is a univariate `fda` object:
  stopifnot(sapply(fd_obj_list, function(x) {
    length(dim(x$coefs)) == 2
  }))
  

  # -------------------------------------------------------------------------
  # extract univariate `fda` objects and elements
  fd_obj_1 <- fd_obj_list[[1]]
  fd_obj_2 <- fd_obj_list[[2]]
  
  basis_obj_1 <- fd_obj_1$basis
  coef_mat_1 <- fd_obj_1$coef
  
  basis_obj_2 <- fd_obj_2$basis
  coef_mat_2 <- fd_obj_2$coef
  
  K_1 <- nrow(coef_mat_1)
  K_2 <- nrow(coef_mat_2)
  
  # construct concatenated coefficient matrix, combining the two dimensions:
  coef_mat_c <- rbind(coef_mat_1, coef_mat_2)
  
  # Checking number of observations
  if(ncol(coef_mat_1) != ncol(coef_mat_2)) stop("pair of univariate functional observations must have the same number of observations")
  N <- ncol(coef_mat_1)
  if (N < 2) stop("PCA not possible without replications.")
  # Check the id_vec:
  stopifnot(is.integer(id_vec))
  if(length(id_vec) != N) stop("Length of id_vec must match number of observations in the fd_obj")
  
  N_subject <- length(unique(id_vec))
  if(!(N_subject < N)) stop("Must have replicate observations from at least some subjects")
  
  # Create design matrix:
  id_factor <- factor(id_vec)
  Zmat <- sparse.model.matrix(~ - 1 + id_factor)
  stopifnot(dim(Zmat) == c(N, N_subject))
  
  # And put into list with indetity matrix for our covariance estimation function (requires this because based of denseFLMM)
  Z_list <- list(list(), list())
  Z_list[[1]][[1]] <- Zmat
  Z_list[[2]][[1]] <- Diagonal(n = N)
  
  # Do a least squares fit on the coefficients (Cederbaum 2017 Section 3/ Appendix A of our paper)
  K_coef_list <- cov_unstruct_coef(Ycoef = t(coef_mat_c),
                                   Zlist = Z_list, 
                                   G = 1, 
                                   Lvec = c(N_subject, N),
                                   groups = NA,
                                   Zvars = NA)
  
  stopifnot(length(K_coef_list) == 2)
  
  # covariance functions are now matrix valued
  K_U_coef <- K_coef_list[[1]] # get basis coefs for the between-subject covariance function
  
  stopifnot(dim(K_U_coef) == rep(K_1 + K_2, 2)) # get basis coefs for the between-subject covariance function

  # # let's make 
  K_U_bifd_list <- list()
  K_U_bifd_list[["K_U_11"]] <- bifd(coef = K_U_coef[(1:K_1), (1:K_1)], # and use coefficients to create bivariate fd object
                                     sbasisobj = basis_obj_1,
                                     tbasisobj = basis_obj_1)
  
  K_U_bifd_list[["K_U_12"]] <- bifd(coef = K_U_coef[(1:K_1), ((K_1 + 1):(K_1 + K_2))], # and use coefficients to create bivariate fd object
                                    sbasisobj = basis_obj_1,
                                    tbasisobj = basis_obj_2)
  K_U_bifd_list[["K_U_21"]] <- bifd(coef = K_U_coef[((K_1 + 1):(K_1 + K_2)),(1:K_1)], # and use coefficients to create bivariate fd object
                                    sbasisobj = basis_obj_2,
                                    tbasisobj = basis_obj_1)
  
  K_U_bifd_list[["K_U_22"]] <- bifd(coef = K_U_coef[((K_1 + 1):(K_1 + K_2)),((K_1 + 1):(K_1 + K_2))], # and use coefficients to create bivariate fd object
                                    sbasisobj = basis_obj_2,
                                    tbasisobj = basis_obj_2)
  # 
  # 
  # 
  # 
  K_E_coef <- K_coef_list[[2]] # get basis coefs for the within-subject covariance function
  stopifnot(dim(K_E_coef) == rep(K_1 + K_2, 2))

  K_E_bifd_list <- list()
  K_E_bifd_list[["K_E_11"]] <- bifd(coef = K_E_coef[(1:K_1), (1:K_1)], # and use coefficients to create bivariate fd object
                                    sbasisobj = basis_obj_1,
                                    tbasisobj = basis_obj_1)
  
  K_E_bifd_list[["K_E_12"]] <- bifd(coef = K_E_coef[(1:K_1), ((K_1 + 1):(K_1 + K_2))], # and use coefficients to create bivariate fd object
                                    sbasisobj = basis_obj_1,
                                    tbasisobj = basis_obj_2)
  K_E_bifd_list[["K_E_21"]] <- bifd(coef = K_E_coef[((K_1 + 1):(K_1 + K_2)),(1:K_1)], # and use coefficients to create bivariate fd object
                                    sbasisobj = basis_obj_2,
                                    tbasisobj = basis_obj_1)
  
  K_E_bifd_list[["K_E_22"]] <- bifd(coef = K_E_coef[((K_1 + 1):(K_1 + K_2)),((K_1 + 1):(K_1 + K_2))], # and use coefficients to create bivariate fd object
                                    sbasisobj = basis_obj_2,
                                    tbasisobj = basis_obj_2)
  # Return list:
  list(K_U_bifd_list = K_U_bifd_list,
       K_E_bifd_list = K_E_bifd_list)
}