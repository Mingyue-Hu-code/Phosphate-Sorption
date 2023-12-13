gauss.newton<-function(xmat, y, ps, fctn, ders, wts)
  {
    # fctn is a function that computes the current 
    # expectations, ders is a function that computes
    # the current matrix of derivatives of the 
    # expectation function (n by p matrix), wts is a function
    # that computes weights as 1/g^2 where g is
    # the variance model
    # xmat are covariates that can be in any
    # form (vector, matrix, list, etc.) that is
    # expected by fctn, ders, and wts
    cps <- ps
    cnt <- 0
    crit <- 1e-08
    repeat {
      cnt <- cnt + 1
      cps1 <- cps
      ee <- fctn(xmat, cps)
      V <- ders(xmat, cps)
      #print(V)
      W <- wts(xmat, cps)
      #print(W)
      d <- solve(t(V) %*% W %*% V)
      #cat("d:",d,fill=T)
      d <- d %*% t(V) %*% W %*% as.matrix(y - ee)
      cps <- cps1 + d
      cat("New Estimates at Iteration",cnt,":", fill = T)
      cat(cps, fill = T)
      dist <- (sum((cps - cps1)^2))^0.5
      if(dist < crit) break
    }
    cat("Convergence Criterion of ", crit, " met", fill = T)
    cat("Final Estimates: ", cps, fill = T)
    cps
  }

nonlin<-function(xmat, ys, ps, fctn, ders, wts)
  {
    # external function definitions for fctn,
    # ders and wts are as in gauss.newton
    # output is list containing betahat(bs),
    # sigma squared hat (sshat), the covariance
    # matrix for betahat (cov), fitted values
    # (yh), studentized residuals (bi), 
    # absolute studentized residuals to
    # (2/3) power (abi),
    # and the matrix of derivatives of the
    # expectation function (fb)
    N <- length(ys)
    P <- length(ps)
    bs <- gauss.newton(xmat, ys, ps, fctn, ders, wts)
    yh <- fctn(xmat, bs)
    r <- ys - yh
    w <- wts(xmat, bs)
    g2 <- 1/(diag(w))
    ri <- r/(g2^0.5)
    fb <- ders(xmat, bs)
    G <- matrix(sqrt(g2), N, P, byrow = F)
    xx <- fb/G
    H <- xx %*% (solve(t(xx) %*% xx)) %*% t(xx)
    h <- diag(H)
    sshat <- (1/(N - P)) * sum(ri^2)
    cov <- matrix(0, P, P)
    cnt <- 0
    repeat {
      cnt <- cnt + 1
      tfb <- as.matrix(fb[cnt,  ])
      tfbfb <- tfb %*% t(tfb)
      tel <- tfbfb/g2[cnt]
      cov <- cov + tel
      if(cnt == N)
        break
    }
    cov <- sshat * solve(cov)
    bi <- ri/((sshat * (1 - h))^0.5)
    abi <- (abs(bi))^(2/3)
    result <- list(bs=bs, sshat=sshat, covb=cov, yhat=yh, stdres=bi, absres=abi, fb=fb)
    result
  }
