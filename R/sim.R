#' Simulation of the mfBm
#'
#' @description
#' This function simulates a multivariate fractional Brownian motion of dimension \eqn{d}
#' using circulant embedding. The input parameters must satisfy the condition
#' for the existence of the covariance matrix, which is verified by the function.
#'
#' @author Prof. J.F. Coeurjolly
#'
#' @param n Length of the time series.
#' @param H Vector of Hurst coefficients.
#' @param sig Vector of scale coefficients.
#' @param rho Symmetric matrix of instantaneous correlation coefficients.
#' @param eta Antisymmetric matrix of asymmetry coefficients.
#' @param delta Time step between consecutive observations in the time series y (e.g., 1/252 for daily data if time is in years).
#' @param m Number of simulations.
#' @param gn Logical. If TRUE, the increments process is returned (multivariate fractional Gaussian noise).
#'
#' @returns list of m numeric matrices of size  d x n
#' @export
#'
#' @references
#' Amblard, Pierre-Olivier, et al. "Basic properties of the multivariate fractional Brownian motion." arXiv preprint arXiv:1007.0828 (2010).
sim_mfbm <- function(n=50, H=c(.2,.2),
                     sig = c(1, 1),
                     rho = matrix(c(1, .5, .5, 1), nc = 2),
                     eta = matrix(c(0, .5, -.5, 0), nc = 2),
                     delta = 1, m = 1, gn = FALSE,
                     plot = FALSE, print = TRUE, choix = NULL, forceEta = FALSE){

  ## Simulation of a multivariate fractional Brownian motion.
  ## m : number of sample paths if =1 a (n,p) matrix is returned
  ## if >1 a list of (n,p) matrices is returned
  ##  notation: causal, wb, eta stand for the causal MFBM, the well-balanced MFBM,
  ## a general MFBM with prescribed eta matrix.
  ## forceEta=TRUE means that eta_jk is changed by eta_jk/(1-H[j]-H[k]) only when H[j]+H[k] !=1
  ## so that this is easy to check the semidefinite positiveness using ellipsis in
  ## See the reference for details.
  ##
  ## Reference. Amblard, Coeurjolly, Philippe and Lavancier (2012)
  ## Basic properties of the multivariate fractional Brownian motion.
  ## Bulletin de la SMF, SC)minaires et CongrC)s, 28, 65-87.
  ns = m
  p<-length(H)
  G<-matrix(0,nr=p,nc=p)
  if (is.null(choix)) choix<-"eta"
  switch(choix,
         "eta"={
           if (forceEta){
             for (j in (1:p)) { for (k in (1:p)){
               if((H[j]+H[k])!=1 ) eta[j,k]<-eta[j,k]/(1-H[j]-H[k])

             }}
           }
         },
         "causal"={
           eta<-matrix(0,nc=p,nr=p)
           for (j in (1:p))
             for (k in (1:p)){
               Hjk<-H[j]+H[k]
               if (Hjk!=1){ eta[j,k]<- rho[j,k]*(cos(pi*H[j])-cos(pi*H[k]))/(cos(pi*H[j])+cos(pi*H[k]))	}
               else{ eta[j,k]<-rho[j,k]*2/pi/tan(pi*H[j])}
             }
           diag(eta)<-0
         },
         "wb"={
           eta<-matrix(0,nc=p,nr=p)
         })

  ## computation of the matrix G (existence proposition)
  for (j in (1:p)){for (k in (1:p)){
    Hjk<-H[j]+H[k]
    if (Hjk!=1){
      pre<-rho[j,k]*sin(pi/2*Hjk)
      pim<- -eta[j,k]*cos(pi/2*Hjk)
    }
    else{
      pre<-rho[j,k];pim<- -pi/2*eta[j,k]
    }
    G[j,k]<-complex(real=pre,imaginary=pim)
    G[j,k]<-gamma(H[j]+H[k]+1)*G[j,k]
  }}
  ltmp<-sum(eigen(G)$values<0)
  #print(H)
  #ltmp<-existMFBM(H=H,rho=rho,eta=eta,forceEta=forceEta,choix=choix)
  txt<-paste("cannot be simulated for these parameters :",ltmp," negative eigenvalues")
  if (print & (ltmp==0)) cat('No negative eigenvalues for G \n')
  if (ltmp>0) stop(txt)

  xlogx<-function(h){
    res<-NULL
    for (e in h) { if (e==0) res<-c(res,0) else res<-c(res,h*log(abs(h)))}
    res
  }

  gammauv<-function(h,Hj,Hk,rhojk,etajk){
    #pour h>=1
    Hjk<-Hj+Hk
    wjk<-function(h,rhojk,etajk){
      tmp2<-NULL
      for (e in h){
        if (Hjk!=1) {
          tmp<-(rhojk-etajk*sign(e))*abs(e)^Hjk
        }	else { tmp<-rhojk*abs(e)-etajk*xlogx(e)}
        tmp2<-c(tmp2,tmp)
      }
      tmp2
    }
    .5*(wjk(h-1,rhojk,etajk)-2*wjk(h,rhojk,etajk)+wjk(h+1,rhojk,etajk))
  }
  cuvj <- function(u,v, m){
    ## provides the first line of the matrix C_m(H_u,H_v) de Wood et Chan version p>1 !
    Hu<-H[u]
    Hv<-H[v]
    Huv<-Hu+Hv
    z0<-rho[u,v]
    z1<-gammauv(1:(m/2-1),Hu,Hv,rho[u,v],eta[u,v])
    z2<-gammauv(m/2,Hu,Hv,rho[u,v],eta[u,v])+gammauv(m/2,Hv,Hu,rho[v,u],eta[v,u])
    z2<-z2/2
    z3<-gammauv(m-(m/2+1):(m-1),Hv,Hu,rho[v,u],eta[v,u])
    c(z0,z1,z2,z3)
  }
  m<-2^(trunc(log(n)/log(2))+2)
  tab<-B<-array(0,dim=c(p,p,m))
  for (u in 1:p){ for (v in 1:p){
    vpCuv <- cuvj( u, v,m)
    vpCuv <- (fft(c(vpCuv), inverse = FALSE))
    vpCuv<-(vpCuv)
    tab[u,v,]<-vpCuv
  }}
  res<-list()
  for (iSamp in (1:ns)){
    cat('\r i=',iSamp)
    W<-matrix(0,nr=m,ncol=p)
    ## Simulations of U and V and storage in an array
    Z<-matrix(0,nr=m,nc=p)
    Z[1,]<-rnorm(p)/sqrt(m)
    Z[m/2+1,]<-rnorm(p)/sqrt(m)

    trace<-NULL
    for (j in (0:(m-1))){
      if ((j>0) & (j<=(m/2-1))){
        U<-rnorm(p); V<-rnorm(p)
        Z[j+1,]<-complex(real=U,imaginary=V) /sqrt(2*m)
        Z[m-j+1,]<-complex(real=U,imaginary=(-V)) /sqrt(2*m)
      }
      A<-tab[,,j+1]
      tmp<-eigen(A)
      vpA<-Re(tmp$values)
      vecpA<-	tmp$vec
      vecpAinv<-solve((tmp$vec))
      vpAPlus<-apply(cbind(vpA,rep(0,length(vpA))),1,max)
      vpANeg<-apply(cbind(-vpA,rep(0,length(vpA))),1,max)
      trace<-c(trace,sum(vpANeg)/sum(vpA))
      B[,,j+1]<-vecpA %*% diag(sqrt((vpAPlus))) %*% vecpAinv
      W[j+1,]<-B[,,j+1] %*% Z[j+1,]
    }
    X <- t(t(Re(mvfft(W,inverse=F)[1:n,])) * (sig * delta  ^ H))

    if(gn == FALSE) {
      vFBM<-apply(X,2,cumsum)
      res[[iSamp]]<-vFBM
    } else if (gn == TRUE) {
      res[[iSamp]]<-X
    }
    # if (print) cat('  Mean trace:',mean(trace),'\n')
  }
  if (ns==1) {  res2<- res[[1]]} else res2<-res
  res2

}

#' Simulation of the mfOU process
#'
#' @description
#' These functions simulate the multivariate fractional Ornstein-Uhlenbeck process of
#' dimension \eqn{d}.
#'
#' - `sim_mfou_app` relies on a combination of circulant embedding and the Euler-Maruyama scheme,
#' and due to the latter step is approximate;
#' - `sim_mfou_exa` uses the Cholesky decomposition of the covariance matrix,
#' and is therefore exact. This methodology might lead to numerical instability
#' for long series or certain parameters' combinations.
#'
#' The input parameters \eqn{\rho}, \eqn{\eta}, and \eqn{H} must satisfy the condition
#' for the existence of the covariance matrix of the underlying mfBm, which is verified
#' at the beginning.
#'
#' @param n Length of the time series.
#' @param H Vector of Hurst coefficients.
#' @param a Vector of speed of mean reversion coefficients.
#' @param nu Vector of diffusion coefficients.
#' @param mu Vector of long-term means.
#' @param rho Symmetric matrix of instantaneous correlation coefficients.
#' @param eta Antisymmetric matrix of asymmetry coefficients.
#' @param delta Time step between consecutive observations in the time series y (e.g., 1/252 for daily data if time is in years).
#' @param m Number of simulations.
#'
#' @returns List of m numeric matrices of size n x d.
#' @export
#' @rdname sim_mfou
#'
#' @examples
#' # parameters
#' m = 3; delta = 1/252; n = 4 / delta
#' H = c(0.1, 0.2, 0.3)
#' rho = matrix(c(1, .3, .4, 0.3, 1, 0.5, 0.4, 0.5, 1), nc = 3)
#' eta = matrix(c(0, .02, 0, - 0.02, 0, - 0.01, 0, 0.01, 0), nc = 3)
#' a = c(1, 2, 3); nu = c(1, 2, 3); mu = c(0, 0, 0)
#'
#' # we simulate lambda observations for each step of size delta
#' # so to reduce the discretization error of the Euler-Maruyama scheme
#' lambda = 24
#'
#' # simulation
#' tmp = sim_mfou_app(n * lambda, delta / lambda , H, rho, eta, a, nu, mu, m)
#'
#' # retaining observations at a distance delta and discarding
#' # the first n/4 observations to approach the stationary distribution
#' mfou = lapply(X = tmp, FUN = function(x) x[seq(lambda, nrow(x), by = lambda), ][(n / 4 + 1) : n, ])
#'
#' # exact simulation
#' sim_mfou_exa(n, H, rho, eta, a, nu, mu, delta)
sim_mfou_app <- function(n, H, rho, eta, a, nu, mu,
                     delta = 1, m = 1) {
  d = length(a)
  mfgn = sim_mfbm(n - 1, H, rep(1, d),
                  rho, eta, gn = TRUE, m = m)
  euler_iter <- function(mfgn){
    mfou = matrix(NA, nrow = n, ncol = d)
    mfou[1, ] = mapply(FUN = rnorm, n = 1, mean = mu, sd = sqrt(var_fou(a, nu, H)))
    for(t in 2 : n) {
      mfou[t, ] = mfou[t - 1, ] + a * (mu - mfou[t - 1, ]) * delta +
        nu * delta ^ H * mfgn[t - 1, ]
    }
    return(mfou)
  }
  if(m == 1){
    mfou = euler_iter(mfgn)
  } else if (m > 1) {
    mfou = lapply(X = mfgn, FUN = euler_iter)
  }
  return(mfou)
}

#' @export
#' @rdname sim_mfou
sim_mfou_exa <- function(n, H, rho, eta, a, nu, mu, delta, m = 1) {
  if(any(cov_existence_mfbm(H, rho, eta, out = TRUE) < 0 )) stop()
  d = length(H)
  if(!all(nrow(rho) == ncol(eta), ncol(rho) == nrow(eta),
          ncol(eta) == d)) stop("Error: the parameters dimensions doesn't match.")

  syy = cov_dt_mfou(n = n, a = a, nu = nu, H = H, rho = rho, eta = eta, delta = delta)
  syy_lt = t(chol(syy))
  tmp = syy_lt %*% matrix(rnorm(n = n * d * m), nr = n * d, nc = m) + mu
  if(m == 1) {
    mfou = tmp[seq(1, length(tmp), d), ]
    for(i in 2 : d) {
      mfou = cbind(mfou, tmp[seq(i, length(tmp), d), ])
    }
    colnames(mfou) <- paste0("y", 1 : d)
  } else if (m > 1) {
    mfou = list()
    for(i in 1 : d) {
      mfou[[i]] <- tmp[seq(i, nrow(tmp), d), ]
    }
    names(mfou) <- paste0("y", 1 : d)
  }
  return(mfou)
}

#' @export
#' @rdname sim_mfou
sim_fou_exa <- function(n, H, a, nu, mu, delta, m = 1) {
  syy = cov_dt_fou(n = n, a = a, nu = nu, H = H, delta = delta)
  syy_lt = t(chol(syy))
  tmp = syy_lt %*% matrix(rnorm(n = n * m), nr = n, nc = m) + mu
  return(tmp)
}
