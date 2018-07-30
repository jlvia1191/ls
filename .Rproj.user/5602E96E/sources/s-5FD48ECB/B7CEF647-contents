# ls.R

#' @title Estimation of the log Likelihood of the Saturated Model
#' @description When the values of the outcome variable Y are either 0 or 1, the function lsm() calculates the estimation of the log likelihood in the saturated model. This model is characterized by Llinas (2006, ISSN:2389-8976) in section 2.3 through the assumptions 1 and 2. If Y is dichotomous and the data are grouped in J populations, it is recommended to use the function lsm() because it works very well for all K.

#' @param formula An expression of the form y ~ model, where y is the outcome variable (binary or dichotomous: its values are 0 or 1).
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lsm() is called.
#' @return  \code{ls} returns an object of class  \code{"ls"}.
#'
#'An object of class  \code{"ls"} is a list containing at least the following components:
#'
#' \item{log_Likelihood}{Estimation of the log likelihood.}
#' \item{populations}{Total number J of populations in the model.}
#' \item{z_j}{ Value of Zj (the sum of the observations in the jth population)}.

#' @details The saturated model is characterized by the assumptions 1 and 2 presented in section 2.3 by Llinas (2006, ISSN:2389-8976).
#' @references [1] Humberto Jesus Llinas. (2006). Accuracies in the theory of the logistic models. Revista Colombiana De Estadistica,29(2), 242-244.
#' @references [2] Hosmer, D. (2013). Wiley Series in Probability and Statistics Ser. : Applied Logistic Regression (3). New York: John Wiley & Sons, Incorporated.
#' @author Humberto Llinas Solano [aut], Universidad del Norte, Barranquilla-Colombia \\ Omar Fabregas Cera [aut], Universidad del Norte, Barranquilla-Colombia \\ Jorge Villalba Acevedo [cre, aut], Unicolombo, Cartagena-Colombia.

#' @export
#' @import stats

ls <- function(formula , data )
{
  mf <- model.frame(formula = formula, data = data)

  res <-do.call(rbind, (tapply(as.vector(mf[, 1]), t(apply((mf[, -1,drop =FALSE]), 1, paste0,collapse = "-")),function(x) c(z = sum(as.numeric(x)), n = length(as.numeric(x)),p = mean(as.numeric(x))))))
  zj<- res[, 1]; nj <- res[, 2]; pj <- res[, 3]; vj <- pj*(1-pj); mj <- nj*pj; Vj <- nj*vj; V <- diag(vj);sp <- as.matrix((zj - nj * pj)/vj); ip <- diag(nj/vj); Zj <- (zj - nj*pj)/sqrt(nj*vj)
  sj <- (res[, 1] * log(res[, 3]) + (res[, 2] - res[, 1]) * log(1 - res[, 3]))
  Lj <-ifelse((res[, 3]) == 0 | (res[, 3]) == 1, 0, sj)
  sat <- sum (Lj)
r<-list(log_Likelihood = sat, populations = length(res) / 3,z_j = as.matrix(zj), n_j = nj, p_j = pj, fitted.values = Lj, v_j = vj, m_j = as.matrix(mj), V_j = Vj, V = V, S_p = sp, I_p = ip, Zast_j = as.matrix(Zj))


}

lsm<- function(x, ...) UseMethod("lsm")

lsm.default <- function(formula , data)
{
  est <- ls(formula , data)
  est$call <- match.call()
  class(est) <- "ls"
}

print.lsm  <- function(x, ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\nLog_Likelihood: \n")
  print(x$log_Likelihood)
  cat("\nPopulations: \n")
  print(x$populations)
}







