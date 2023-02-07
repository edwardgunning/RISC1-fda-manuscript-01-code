project_mean_onto_fpcs <- function(pca.fd_obj) {
  # adapted code from `fda` package for projecting
  # individual functions #onto FPCs to project mean
  # onto FPCs.
  # The solo argument pca.fd_obj must be a
  # (uni- or multivariate) pca.fd object.
  stopifnot(class(pca.fd_obj) == "pca.fd")
  
  if(sum(pca.fd_obj$varprop) <= 0.999) {
    warning("pve <= 0.99: projecting the mean onto the FPCs may not be a good idea if transformation is not (near-)lossless!")
    }
  harmfd <- pca.fd_obj$harmonics
  stopifnot(is.fd(harmfd))
  
  mean_fd_obj <- pca.fd_obj$meanfd
  stopifnot(is.fd(mean_fd_obj))

  mean_coefs <- mean_fd_obj$coefs
  
  ndim  <- length(dim(mean_coefs))
  nharm <- dim(harmfd$coefs)[2]
  
  if(ndim == 2) {
    nvar <- 1
  } else if(ndim == 3) {
    nvar <- dim(mean_coefs)[3]
  } else stop("Dimension of coefficient matrix wrong!")
  
  
  if (nvar == 1) {
    harmscr  <- inprod(mean_fd_obj, harmfd)
  } else {
    harmscr  <- array(NA, c(1,   nharm,  nvar))
    coefarray <- mean_fd_obj$coefs
    harmcoefarray <- harmfd$coefs
    for (j in 1:nvar) {
      fdobjj  <- fd(as.matrix(coefarray[,,j]), mean_fd_obj$basis)
      harmfdj <- fd(as.matrix(harmcoefarray[,,j]), harmfd$basis)
      harmscr[,,j] <- inprod(fdobjj, harmfdj)
    }
  }
  
  harmscr
}



