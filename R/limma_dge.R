limma_dge <- function(expression_data,caseIds=list(),ctrIds) {
  
  if(!is.list(caseIds)) #ie, only one case
  {
    tmp_data <- expression_data[,c(ctrIds,caseIds)]
    
    expdes <- matrix(0,nrow=ncol(tmp_data),2) # experimental design -- initialize matrix with 0s -- these will eventually become one-hot encodings for case vs. control
    expdes[seq(1,length(ctrIds)),1] <- 1
    expdes[seq(length(ctrIds)+1,length(ctrIds)+length(caseIds)),2] <- 1
    colnames(expdes) <- c("control","case")
    
    contmat <- makeContrasts(case-control,levels=expdes) # contrast matrix
    
    limma.fit <- lmFit(tmp_data,design=expdes)
    limma.fit <- contrasts.fit(limma.fit,contmat)
    limma.fit <- eBayes(limma.fit)
    
    res <- topTable(limma.fit,1,number=nrow(tmp_data),sort.by = "none",adjust.method = "fdr")
    return(res)
    
  } else {
    # one-hot encode
    expdes <- matrix(0,nrow=ncol(expression_data),1) 
    expdes[ctrIds,1] <- 1 # controls
    
    for(i in seq(1,length(caseIds)))
    {
      x <- matrix(0,nrow=ncol(expression_data),1) 
      x[caseIds[[i]],1] <- 1
      expdes <- cbind(expdes,x)
    }
    
    tmp_data <- expression_data[,c(ctrIds,unlist(caseIds))]
    expdes <- expdes[c(ctrIds,unlist(caseIds)),]
    
    colnames(expdes) <- c("control",paste("case",seq(1,length(caseIds)),sep=""))
    myContrasts <- c(paste(colnames(expdes)[2:ncol(expdes)],"-",colnames(expdes)[1],sep=""))
    contmat <- contmat <- eval(as.call(c(as.symbol("makeContrasts"),as.list
                                         (myContrasts),levels=list(expdes)))) # contrast matrix
    
    limma.fit <- lmFit(tmp_data,design=expdes)
    limma.fit <- contrasts.fit(limma.fit,contmat)
    limma.fit <- eBayes(limma.fit)
    
    return(limma.fit)
  }
}