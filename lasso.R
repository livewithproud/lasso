# diabetes <- read.table("/home/scuiting/Documents/diabetes.txt", header=TRUE)
# X <- data.matrix(diabetes[,1:10])
# Y <- data.matrix(diabetes[,11])


lasso <- function(x , y, normalize = TRUE, intercept = TRUE, eps = .Machine$double.eps, max.steps, trace = TRUE)
{
  call <- match.call()
  ### 参数
  nm <- dim(x)
  n <- nm[1]
  m <- nm[2]
  one <- rep(1,n)
  im <- inactive <- seq(m)
  vn <- dimnames(x)[[2]]

  ### x的中心化，标准化   y 的中心化
  if(intercept){
    meanx <- drop(one %*% x)/n
    x <- scale(x, meanx, FALSE)
    mu <- mean(y)
    y <- drop(y-mu)
  }else {
    meanx <- rep(0,m)
    mu <- 0
    y <- drop(y)
  }
  
  if(normalize){
    normx <- sqrt(drop(one %*% (x^2)))
    # 无效变量
    nosignal<-normx/sqrt(n) < eps
    if(any(nosignal))# ignore variables with too small a variance
    {
      ignores<-im[nosignal]
      inactive<-im[-ignores]
      normx[nosignal]<-eps*sqrt(n)
      if(trace){
        cat("LARS Step 0 :\t", sum(nosignal), "Variables with Variance < \ eps; dropped for good\n")  #
      }
    }else ignores <- NULL
  
    names(normx) <- NULL
    x <- scale(x, FALSE, normx)
  }else {
    normx <- rep(1,m)
    ignores <- NULL
  }  
 
  
  ### 相关系数向量
  ssy <- sum(y^2)
  Cvec <- drop(t(y) %*% x)  
  # 参数
  residuals <- y
  if(missing(max.steps)) max.steps <- 8*min(m, n-intercept)
  beta <- matrix(0, max.steps+1, m)
  lambda <- double(max.steps)
  R2 <- 1
  RSS <- ssy
  first.in <- integer(m) # m个0
  active <- NULL
  actions <- as.list(seq(max.steps))
  drops <- FALSE
  Sign <- NULL
  R <- NULL
  k <- 0
  ### 函数正文
  while((k < max.steps) & (length(active) < min(m - length(ignores), n - intercept)))
  {
    action <- NULL
    C <- Cvec[inactive]
    Cmax <- max(abs(C))
    if(Cmax<eps*100){ # the 100 is there as a safety net
      if(trace) cat("Max |corr| = 0; exiting...\n")
      break
    }
    k <- k+1
    lambda[k] <- Cmax
    if(!any(drops)){  # 如果drop==FALSE
      new <- abs(C) >= Cmax - eps
      C <- C[!new]
      new <- inactive[new]
      for(inew in new) { # active Sign action
        R <- updateR(x[, inew], R, x[, active])
        if(attr(R, "rank") == length(active)) {
          nR <- seq(length(active)) 
          R <- R[nR, nR, drop = FALSE]
          attr(R, "rank") <- length(active)
          ignores <- c(ignores, inew)
          action <- c(action, -inew)
          if(trace) {
            cat("LASSO Steps", k, ":\t Variable", inew, "\t collinear; dropped for good\n")
          }
        }
        else{
          if(first.in[inew] == 0) first.in[inew] <- k
          active <- c(active, inew)
          Sign <- c(Sign, sign(Cvec[inew]))
          action <- c(action, inew)
          if(trace) {
            cat("\nLASSO Steps", k, ":\t Variable", inew, "\tadded\n")
          }
        }      
      }
    }
    else action <- -dropid
    Gi1 <- backsolve(R, backsolve(R, Sign, ncol(R), transpose=TRUE))
    A <- 1/sqrt(sum(Gi1 * Sign))
    w <- A*Gi1  #  with direction
    u <- drop(x[,active, drop=FALSE] %*% w)
    if(length(active) >= min(n-1, m-length(ignores))) {
      gamhat <- Cmax/A
    }
    else {
      a <- drop(u %*% x[, -c(active, ignores), drop = FALSE])
      gam <- c((Cmax - C)/(A-a), (Cmax + C)/(A+a))
      gamhat <- min(gam[gam > eps], Cmax/A)
    }    
    # 是否减去变量
    dropid <- NULL
    b1 <- beta[k, active]
    z1 <- -b1/w   
    zmin <- min(z1[z1 > eps], gamhat)
    if(zmin < gamhat) {
      gamhat <- zmin
      drops <- z1 == zmin
    }
    else drops <- FALSE
    # 系数
    beta[k + 1, ] <- beta[k, ]
    beta[k + 1, active] <- beta[k + 1, active] + gamhat * w
    # residuals
    residuals <- residuals - gamhat * u   
    Cvec <- drop(t(residuals) %*% x)

    if(any(drops)) {
      dropid <- seq(drops)[drops]
      for(id in rev(dropid)) {
        if(trace) {
          cat("Lasso Step", k+1, ":\t Variable", active[id], "\tdropped\n")
        }
        R <- downdateR(R, id)
      }
      dropid <- active[drops]
      beta[k+1,dropid] <- 0
      active <- active[!drops]
      Sign <- Sign[!drops]
    }
    if(!is.null(vn))  names(action) <- vn[abs(action)]  # 如果列名称不为空 
    actions[[k]] <- action
    inactive <- im[ - c(active,ignores)]
  }
  beta <- beta[seq(k + 1), ,drop = FALSE]
  lambda <- lambda[seq(k)]
  dimnames(beta) <- list(paste(0:k), vn)
  # RSS, R2,Cp
  if(trace) {
  cat("Computing residuals, RSS etc .....\n")
  }
  residuals <- y - x %*% t(beta)
  beta <- scale(beta, FALSE, normx)
  RSS <- apply(residuals^2, 2, sum)
  R2 <- 1 - RSS/RSS[1]
  actions <- actions[seq(k)]
  netdf <- sapply(actions, function(x) sum(sign(x)))
  df <- cumsum(netdf)
  if(intercept) {
    df <- c(Intercept <- 1,df+1)
  }else { df <- c(NULL <- 0, df)}
  rss.big <- rev(RSS)[1]
  df.big <- n-rev(df)[1]
  if(rss.big < eps | df.big < eps) sigma2 <- NaN
  else
    sigma2 <- rss.big/df.big
  Cp <- RSS/sigma2 - n + 2 * df ####### 
  object <- list(beta = beta, Cp = Cp, RSS = RSS, df = df,
                 netdf = netdf,actions=actions,TSS = ssy)
  ### R2 = R2, A = A, ignores <- ignores
  ### actions = actions[seq(k)],entry = first.in,mu = mu, normx = normx, meanx = meanx, x=x, y=y
  object
}

updateR <- function(xnew, R = NULL, xold, eps = .Machine$double.eps)
{
  ###Gram argument determines the nature of xnew and xold
  xtx <- sum(xnew^2)
  norm.xnew <- sqrt(xtx)
  if(is.null(R)) {
    R <- matrix(norm.xnew, 1, 1)
    attr(R, "rank") <- 1
    return(R)
  }
  Xtx <- drop(t(xnew) %*% xold)
  r <- drop(backsolve(R, Xtx, transpose=TRUE))
  rpp <- norm.xnew^2 - sum(r^2)
  rank <- attr(R, "rank")  ### check if R is machine singular
  if(rpp <= eps)
    rpp <- eps
  else {
    rpp <- sqrt(rpp)
    rank <- rank + 1
  }
  R <- cbind(rbind(R, 0), c(r, rpp))
  attr(R, "rank") <- rank
  R
}

downdateR <-function(R, k = p)
{
  p <- dim(R)[1]
  if(p == 1)  return(NULL)
    r <- R[,  - k, drop = FALSE]
    pl <- p-1
  if(k <= pl) 
  {
    for(i in k:pl)
    {

      a <- r[i,i]
      b <- r[i+1,i]
      if(b == 0) break
      if(abs(b) > abs(a))
      {
        tau = a/b
        s = (1 + tau^2)^(-0.5)
        c = s * tau
      }else{
        tau = b/a
        c = (1 + tau^2)^(-0.5)
        s = c * tau
      }
      r[i,i] = c*a + s*b
      r[i+1,i] = s*a - c*b
      for(j in (i+1) : pl)
      {
        if((i+1) > pl) break
        a = r[i,j]
        b = r[i+1,j]
        r[i,j] = c*a + s*b
        r[i+1,j] = s*a - c*b
      }
    }    
  }
    R <- r[-p,,drop = FALSE]
    attr(R ,"rank") <- p-1
    R
}




