#' CIBERSORT R script v1.03
#' Note: Signature matrix construction is not currently available; use java version for full functionality.
#' Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' Requirements:
#'       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#'       install.packages('e1071')
#'       install.pacakges('parallel')
#'       install.packages('preprocessCore')
#'       if preprocessCore is not available in the repositories you have selected, run the following:
#'           source("http://bioconductor.org/biocLite.R")
#'           biocLite("preprocessCore")
#' Windows users using the R GUI may need to Run as Administrator to install or update packages.
#' This script uses 3 parallel processes.  Since Windows does not support forking, this script will run
#' single-threaded in Windows.
#'
#' Usage:
#'       Navigate to directory containing R script
#'
#'   In R:
#'       source('CIBERSORT.R')
#'       results <- CIBERSORT('sig_matrix_file.txt','mixture_file.txt', perm, QN)
#'
#'       Options:
#'       i)  perm = No. permutations; set to >=100 to calculate p-values (default = 0)
#'       ii) QN = Quantile normalization of input mixture (default = TRUE)
#'
#' Input: signature matrix and mixture file, formatted as specified at http://cibersort.stanford.edu/tutorial.php
#' Output: matrix object containing all results and tabular data written to disk 'CIBERSORT-Results.txt'
#' License: http://cibersort.stanford.edu/CIBERSORT_License.txt
#' Core algorithm
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#' @export
CoreAlg <- function(X, y){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}

#' do permutations
#' @param perm Number of permutations
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#' @export
#' 

doPerm <- function(perm, X, Y){
  message(paste0("Running doPerm loop, number of threads used: ",detectCores()))
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist<-list()
  registerDoParallel(round(detectCores()))
  dist<-foreach(itor = 1:perm) %dopar% {  
    #my_func(x)  
  #}
  #while(itor <= perm){
    #print(itor)
    library(limma)
    library(parallel)
    library(ggplot2)
    library(ggpubr)
    library(future.apply)
    library(car)
    library(ridge)
    library(e1071)
    library(preprocessCore)
    library(tcltk)
    library(limma)
    library(parallel)
    library(ggplot2)
    library(ggpubr)
    library(future.apply)
    library(car)
    library(ridge)
    library(e1071)
    library(preprocessCore)
    library(foreach)  
    library(doParallel)
    CoreAlg <- function(X, y){
      
      #try different values of nu
      svn_itor <- 3
      
      res <- function(i){
        if(i==1){nus <- 0.25}
        if(i==2){nus <- 0.5}
        if(i==3){nus <- 0.75}
        model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
        model
      }
      
      if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
        out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
      
      nusvm <- rep(0,svn_itor)
      corrv <- rep(0,svn_itor)
      
      #do cibersort
      t <- 1
      
      while(t <= svn_itor) {
        weights = t(out[[t]]$coefs) %*% out[[t]]$SV
        weights[which(weights<0)]<-0
        w<-weights/sum(weights)
        u <- sweep(X,MARGIN=2,w,'*')
        k <- apply(u, 1, sum)
        nusvm[t] <- sqrt((mean((k - y)^2)))
        corrv[t] <- cor(k, y)
        t <- t + 1
      }
      
      #pick best model
      rmses <- nusvm
      mn <- which.min(rmses)
      model <- out[[mn]]
      
      #get and normalize coefficients
      q <- t(model$coefs) %*% model$SV
      q[which(q<0)]<-0
      w <- (q/sum(q))
      
      mix_rmse <- rmses[mn]
      mix_r <- corrv[mn]
      
      newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
      
    }
    
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)
    
    mix_r <- result$mix_r
    
    #store correlation
    #if(itor == 1) {dist <- mix_r}
    #else {dist <- rbind(dist, mix_r)}
    dist[[itor]]<-mix_r
    return(dist[[itor]])
  }
  stopImplicitCluster()
  dist<-do.call(rbind,dist)
  newList <- list("dist" = dist)
}

#' Main functions
#' @param sig_matrix file path to gene expression from isolated cells
#' @param mixture_file heterogenous mixed expression
#' @param perm Number of permutations
#' @param QN Perform quantile normalization or not (TRUE/FALSE)
#' @export
CIBERSORT <- function(sig_matrix, mixture_file, perm=0, QN=TRUE){
  library(e1071)
  library(parallel)
  library(preprocessCore)
  
  #read in data
  message("Reading expression matrix")
  X <- fread(sig_matrix,header=T,sep="\t",check.names=F,data.table=FALSE)
  rownames(X)=X[,1]
  X=X[,2:ncol(X)]
  Y <- fread(mixture_file, header=T, sep="\t",check.names=F,data.table=FALSE)
  rownames(Y)=Y[,1]
  Y=Y[,2:ncol(Y)]
  X <- na.omit(X)
  Y <- na.omit(Y)
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]
  
  P <- perm #number of permutations
  
  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}
  
  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  
  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  
  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}
  
  #print(nulldist)
  
  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  #print(header)
  #output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999
  
  #iterate through mixtures
  message(paste0("Running CIBERSORT loop, number of threads used: ",detectCores()))
  output=list()
    registerDoParallel(round(detectCores()))
    output<- foreach(itor = 1:mixtures) %dopar% {  
      #foreach(itor = 1:mixtures) %dopar% {  
      #my_func(x)  
      #}
      #while(itor <= perm){
      #print(itor)
      library(limma)
      library(parallel)
      library(ggplot2)
      library(ggpubr)
      library(future.apply)
      library(car)
      library(ridge)
      library(e1071)
      library(preprocessCore)
      library(tcltk)
      library(limma)
      library(parallel)
      library(ggplot2)
      library(ggpubr)
      library(future.apply)
      library(car)
      library(ridge)
      library(e1071)
      library(preprocessCore)
      library(foreach)  
      library(doParallel)
      CoreAlg <- function(X, y){
        
        #try different values of nu
        svn_itor <- 3
        
        res <- function(i){
          if(i==1){nus <- 0.25}
          if(i==2){nus <- 0.5}
          if(i==3){nus <- 0.75}
          model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
          model
        }
        
        if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
          out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)
        
        nusvm <- rep(0,svn_itor)
        corrv <- rep(0,svn_itor)
        
        #do cibersort
        t <- 1
        
        while(t <= svn_itor) {
          weights = t(out[[t]]$coefs) %*% out[[t]]$SV
          weights[which(weights<0)]<-0
          w<-weights/sum(weights)
          u <- sweep(X,MARGIN=2,w,'*')
          k <- apply(u, 1, sum)
          nusvm[t] <- sqrt((mean((k - y)^2)))
          corrv[t] <- cor(k, y)
          t <- t + 1
        }
        
        #pick best model
        rmses <- nusvm
        mn <- which.min(rmses)
        model <- out[[mn]]
        
        #get and normalize coefficients
        q <- t(model$coefs) %*% model$SV
        q[which(q<0)]<-0
        w <- (q/sum(q))
        
        mix_rmse <- rmses[mn]
        mix_r <- corrv[mn]
        
        newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
        
      }
    
    y <- Y[,itor]
    
    #standardize mixture
    y <- (y - mean(y)) / sd(y)
    
    #run SVR core algorithm
    result <- CoreAlg(X, y)
    
    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}
    #pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))
    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    #out <- c(colnames(Y)[itor],w,mix_r,mix_rmse)
    output[[itor]]<-out
    return(output[[itor]])
    #out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    #if(itor == 1) {output <- out}
    #else {output <- rbind(output, out)}
    }
    stopImplicitCluster()
    output<-do.call(cbind,output)
  output<-t(output)
  output<-rbind(header,output)
 # output=as.data.frame(output)
  #colnames(output)=output[1,]
  #rownames(output)=output[,1]
  #output=output[2:nrow(output),2:ncol(output)]
  #View(output)
 # Correlations=output$Correlation
 # print(Correlations)
 # print(nulldist)
  #pva<-list()
 # plan(multisession)
 # future_lapply(1:length(Correlations),function(i){
  #  Correlation=as.numeric(Correlations[i])
  #  print(nulldist[i])
  #  pva[i]=1 - (which.min(abs(nulldist - Correlation)) / length(nulldist))
  #  return(p)
  #  })
  #pvalue=do.call(cbind,pva)
  #print(pva)
  #output<-data.frame(output,"P-value"=pva)
  write.table(output, file="results/imm_infil/CIBERSORT_Results.txt", sep="\t", row.names=F, col.names=F, quote=F)
  message("output file location: results/imm_infil/CIBERSORT_Results.txt")
}

