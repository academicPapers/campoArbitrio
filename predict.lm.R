pred.lm <- function(mod, predict.df, level = .8, 
                    interval = c("confidence", "prediction"), cluster = NULL, 
                    type = c("HC3", "const", "HC", "HC0", "HC1", "HC2", "HC4",
                             "HC4m", "HC5")){
  
  ##Written by Joshua Gubler ~  http://scholar.byu.edu/jgubler
  ##Last updated on 26 June 2014
  ##Updated by Baobao Zhang
  ##This provides an option for robust (including cluster robust) or non-robust 
  ##standard errors
  ##Note: when estimating a polynomial, you must create the quadratic/cubic as a 
  ##separate variable first!!  
  ##This is also the best procedure when estimating logged effects.  However, 
  ##when estimating interaction effects, there is no need to create a separate 
  ##interaction term.
  ##Also note, that for this function to work well, you must input factor 
  ## variables with more than two levels individually (as indiviual dummies).
  ##Updated by Luiz F. P. Droubi on 06/17/2018 (http://droubi.me)
  
  if(missing(predict.df)){ predict1.df <- mod$model }
  require(sandwich, quietly = TRUE)
  require(lmtest, quietly = TRUE)
  tt <- terms(mod)
  Terms <- delete.response(tt)
  m.mat <- model.matrix(Terms,data=(if(missing(predict.df)){predict1.df}else{predict.df}))
  t1 <- mod$model
  fit <- as.vector(m.mat %*% mod$coef)
  
  if (is.null(cluster)) {
    robcov <- vcovHC(mod, type = type)
    se.fit <-  ifelse(interval == "confidence", 
                      sqrt(diag(m.mat%*%robcov%*%t(m.mat))),
                      sqrt(diag(m.mat%*%robcov%*%t(m.mat)) + 1)
    )
  } else { # cluster robust standard errors
    M <- length(unique(cluster))
    N <- length(cluster)
    K <- mod$rank
    # degrees of freedom adjustment
    dfc <- (M/(M-1))*((N-1)/(N-K))
    # meat
    XDX <- crossprod(apply(estfun(mod), 2, function(x) tapply(x, cluster, sum)))/N
    # variance calculation 
    vcovCL <- dfc*sandwich(mod, meat=XDX)
    se.fit <- ifelse(interval == "confidence", 
                     sqrt(diag(m.mat%*%vcovCL%*%t(m.mat))),
                     sqrt(diag(m.mat%*%vcovCL%*%t(m.mat)) + 1)
                     )
  }
  ci.lower <- fit + qnorm((1 - level)/2)*se.fit
  ci.upper <- fit + qnorm(1 - (1 - level)/2)*se.fit
  
  pred.df <- data.frame(predicted.value = fit, 
                        ci.lower = ci.lower, ci.upper = ci.upper)
  nm <-deparse(substitute(mod))
  mdlname1 <- paste(nm,"allpred.df",sep=".")
  mdlname2 <- paste(nm,"pred.df",sep=".")
  if(missing(predict.df)){
    allpred.df <- cbind(m.mat,t1[1],pred.df)
    assign(mdlname1,allpred.df,envir = .GlobalEnv)}
  else{assign(mdlname2,pred.df,envir = .GlobalEnv)}
}
