#' Return error matrices
#'
#' @param mx A error matrix.
#' @param merr Error value for missing data.
#' @return \code{error_mx} returns a list of three matrices.
#' In this matrix,
#' row1 is true_gen 00, row2 is true_gen 01, row3 is true_gen 11.
#' cols 1-3 are obs. genotype (00,01,11) and last col4 is the missing data.
#' 1. 00x01
#' 2. 01x01
#' 3. 11x01
#' @examples
#' #Genotype error matrix
#' gen_error_mat(hom.error=0.02, het.error=0.8)
#'
#' #Genotype error by consideringMendelian segregation rate
#' error_mx(hom.error=0.02, het.error=0.8)
#'
error_probs <- function(mx, merr){

  probs <- vector("list",3)
  ### remember 4th column is missing data NOT SURE THIS SHOULD BE 1!!!
  #AA by AA (1,0,0), Aa (1/2,1/2,0), aa (0,1,0)
  probs[[1]] <- list(cbind(mx*matrix(c(1, 0, 0), nrow = 3, byrow=F, ncol=3), merr),
                     cbind(mx*matrix(c(1/2, 1/2, 0), nrow = 3, byrow=F, ncol=3), merr),
                     cbind(mx*matrix(c(0, 1, 0), nrow = 3, byrow=F, ncol=3), merr) )
  #Aa by AA (1/2,1/2,0), Aa (1/4,1/2,1/4), aa (0,1/2,1/2)
  probs[[2]] <- list(cbind(mx*matrix(c(1/2,1/2,0), nrow = 3, byrow=F, ncol=3),merr),
                     cbind(mx*matrix(c(1/4,1/2,1/4), nrow = 3, byrow=F, ncol=3),merr),
                     cbind( mx*matrix(c(0,1/2,1/2), nrow = 3, byrow=F, ncol=3),merr) )
  #aa by AA (0,1,0), Aa (0,1/2,1/2), aa (1,0,0)
  probs[[3]] <- list(cbind(mx*matrix(c(0,1,0), nrow = 3, byrow=F, ncol=3),merr),
                     cbind(mx*matrix(c(0,1/2,1/2), nrow = 3, byrow=F, ncol=3),merr),
                     cbind(mx*matrix(c(0,0,1), nrow = 3, byrow=F, ncol=3),merr) )

  return(probs)
}

#' @rdname error_probs
gen_error_mat <- function(major.error, het.error, minor.error){
  mx <- matrix(c(1-major.error,major.error/2,major.error/2,het.error/2,
                 1-het.error,het.error/2,minor.error/2,minor.error/2,1-minor.error),
               byrow=T,nrow=3,ncol=3)
  rownames(mx) <- c("g0", "g1", "g2")
  colnames(mx) <- c("ob0", "ob1", "ob2")
  return(mx)
}
