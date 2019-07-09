library(tidyr)
library(dplyr)
library(readxl)
library(tibble)
library(stringr)
library(rafalib)
library(affy)
###################
# R hacks
###################
slicedat <- function(x, thr=0.1) {
  out <- list()
  xl <- quantile(x,thr, na.rm = T)
  xh <- quantile(x, 1-thr, na.rm=T)
  out$id <- c(which(x <= xl), which(x >= xh))
  out$label <- c(rep(0, length(which(x <= xl))), rep(1, length(which(x >= xh))))
  return(out)
}

patfind <- function(pats, mystr){
  # Searches for multiple terms in pats and returns logical f both elements found in vector element
  # OR AND NOT mode accepted. Example OR: a1 a2, AND: b1 b2, NOT: c1 c2
  #   ==> (a1 OR a2) AND (b1 OR b2) but NOT (c1 OR c2) etc
  if(is.character(pats)){
    ix <- Reduce(`&`, lapply(pats, grepl, mystr))
    return(ix)
  }  else if(is.list(pats)){
    #OR<-pats$or; AND<-pats$and; NOT<-pats$not
    ix <- list()
    for(s in c('or','and','not')){
      ix[[s]] <- Reduce('|', lapply(pats[[s]], grepl, mystr))
    }
    ix0 <- ix[['or']]
    if(!is.null(ix[['and']])) ix0 <- ix0 & ix[['and']]
    if(!is.null(ix[['not']])) ix0 <- ix0 & !ix[['not']]
    return(ix0)
  }
}

center_mean <- function(x) {
  ones = rep(1, ncol(x))
  x_mean = rowMeans(x) %*% t(ones)
  x - x_mean
}

center_columns <- function(x) {
  c <- mean(x)
  c + scale(x, scale=F)
}

xfer_dimnames <- function(x,y){
  rownames(x) <- rownames(y)
  colnames(x) <- colnames(y)
  return(x)
}

coefvar <- function(x){
  sd(x)/mean(x)
}

rankcor <- function(x,y,pval=F){
  if(!pval) return(cor(x,y,method='spearman',use='complete.obs'))
  if(pval) {
    temp <- cor.test(x,y,method='spearman',use='complete.obs')
    out <- list(rho=temp$estimate, pval=temp$p.value)
    return(out)
  }
}

topgenes <- function(m, N=500, type='var', out='mat', cov=NULL){
  
  if(type=='var') {
    if(!is.null(cov)) m <- rowMeans(m) + t(resid(lm(t(m) ~ cov)))
    rv <- genefilter::rowVars(m)
  } else if(type=='sum'){
    rv <- rowSums(m)
  }
  rvord <- order(rv, decreasing=T)[1:N]
  
  if(out=='mat'){
    return(m[rvord,])
  } else if(out=='row'){
    return(rvord)
  }
}

coefvarmat <- function(x){
  genefilter::rowSds(x)/rowMeans(x)
}

vprint <- function(x){
  cat(paste0(x, '\n'))
}

pkgTest <- function(y){
  sapply(y, function(x) {
    #http://stackoverflow.com/questions/9341635/check-for-installed-packages-before-running-install-packages
    if (!require(x,character.only = TRUE))
    {
      install.packages(x,dep=TRUE)
      if(!require(x,character.only = TRUE)) stop("Package not found")
    } else {return(NULL)}
  })
}


###################
# PCA tools
###################
longsvd <- function(x, geneload=FALSE){
  y <- center_mean(x)
  s <- svd(y)
  colnames(s$v) <- paste0('pc', c(1:ncol(s$v)))
  rownames(s$v) <- colnames(x)
  if(! geneload) s$u <- NULL
  return(s)
}

eigfrac <- function(x){
  return(x$d^2/sum(x$d^2))
}

plotscree <- function(x, title='Screeplot', ...){
  y <- x$d^2/sum(x$d^2)
  plot(y, xlab='Eigenvalue $#', ylab='Fraction variation explained', main=title, pch=19)
}

# plotpc <- function(x, pc1=1, pc2=2, group=NULL, col='red', title='PC plot', ...){
#   plot(x$v[,pc1], x$v[,pc2], xlab=paste0('PC', pc1), ylab=paste0('PC', pc2), pch=19, main=title, ...)
#   if(!is.null(group)) points(x$v[group, pc1], x$v[group, pc2], col=col, pch=19)
# }

plotpc <- function(x, pc1=1, pc2=2, group=NULL, title='PCA plot', label=F){
  xax <- paste0('pc', pc1)
  yax <- paste0('pc', pc2)
  varex <- round(x$d^2*100/sum(x$d^2))
  xl <- glue::glue("{xax}: {varex[pc1]}% variance")
  yl <- glue::glue("{yax}: {varex[pc2]}% variance")
  out <-  ggplot2::ggplot(data.frame(x$v), aes_string(xax, yax)) + theme_bw() + xlab(xl) + ylab(yl) + ggtitle(title)
  if(!is.null(group)){
    out <- out + geom_point(aes(col=group))
  } else{
    out <- out + geom_point()
  }
  if(label) out <- out+ggrepel::geom_text_repel(aes(label=rownames(x$v)))
  print(out)
  #plot(x$v[,pc1], x$v[,pc2], xlab=paste0('PC', pc1), ylab=paste0('PC', pc2), pch=19, main=title, ...)
  #if(!is.null(group)) points(x$v[group, pc1], x$v[group, pc2], col=col, pch=19)
}
getOutlier <- function(x, method='median', thr=4){
  if(method=='median'){
    mx <- median(x)
    devx <- mad(x)
  } else if(method=='mean'){
    mx <- mean(x)
    devx <- sd(x)
  }
  
  lx <- mx - devx * thr
  ux <- mx + devx * thr
  
  outlier <- which(x < lx | x > ux)
  return(outlier)
}
###################
# Affy tools
###################
affy2sym <- function(x, PM=F, allprobes = T){
  #require(hgu133plus2.db)
  #require(plyr)
  if(allprobes) {
    y <- toggleProbes(hgu133plus2SYMBOL, "all")
  } else {
    y <- hgu133plus2SYMBOL
  }
  # Convert to a list
  #print(y[x])
  yy <- AnnotationDbi::as.list(y[x])
  yy2 <- sapply(yy, function(z) if (length(z) > 1) { paste0(z, collapse = " /// ")} else (z))
  out <- ldply(yy2)
  if(length(out)>0) colnames(out) <- c("Probe Set ID", "Gene Symbol")
  return(out)
}

affy2uniprot <- function(x, PM=F, allprobes = T){
  #require(hgu133plus2.db)
  #require(plyr)
  if(allprobes) {
    y <- toggleProbes(hgu133plus2UNIPROT, "all")
  } else {
    y <- hgu133plus2UNIPROT
  }
  # Convert to a list
  yy <- AnnotationDbi::as.list(y[x])
  yy2 <- sapply(yy, function(z) if (length(z) > 1) { paste0(z, collapse = " /// ")} else (z))
  out <- ldply(yy2)
  if(length(out)>0) colnames(out) <- c("Probe Set ID", "UniprotID")
  return(out)
}

affy2Entrez <- function(x, PM=F, allprobes = T){
  #require(hgu133plus2.db)
  #require(plyr)
  if(allprobes) {
    y <- toggleProbes(hgu133plus2ENTREZID, "all")
  } else {
    y <- hgu133plus2ENTREZID
  }
  # Convert to a list
  yy <- AnnotationDbi::as.list(y[x])
  yy2 <- sapply(yy, function(z) if (length(z) > 1) { paste0(z, collapse = " /// ")} else (z))
  out <- ldply(yy2)
  colnames(out) <- c("Probe Set ID", "Entrez ID")
  return(out)
}

getAffySymbols <- function(){
  y <- toggleProbes(hgu133plus2.db::hgu133plus2SYMBOL, "all")
  return(unlist(as.list(y)) %>% unique)
}

sym2affy <- function(x, PM=F, allprobes = T){
  #require(hgu133plus2.db)
  #require(plyr)
  #require(reshape2)
  if(allprobes) {
    y <- toggleProbes(hgu133plus2ALIAS2PROBE, "all")
  } else {
    y <- hgu133plus2ALIAS2PROBE
  }
  # Convert to a list
  yy <- AnnotationDbi::as.list(y[x]) 
  out <- melt(yy)
  colnames(out) <- c("Probe Set ID", "Gene Symbol")
  return(out)
}

catPaste0 <- function(x) cat(paste0(x,'\n'))

hkg2affy <- function(x){
  hkgenes <- read.csv('~/Documents/work/data/affy/HKG_Hendrik_2007.csv', header=T, as.is=T)
  pr <- sym2affy(hkgenes$Genesymbol)
  return(levels(pr$`Probe Set ID`)[unclass(pr$`Probe Set ID`)])
}

affy_gender_caller <- function(datn, PM=F){
  #Gender caller
  tag <- ifelse(PM, '_PM', '')
  mp<- paste0('201909', tag, '_at')
  fp <- paste0('224588', tag, '_at')
  gender_dat <- datn[c(mp, fp), ]
  gender <- apply(gender_dat, 2, function(x) {if(x[2]/x[1]>2){"Female"}else if((x[2]/x[1]>0.7)&(x[2]/x[1]<=2)){"Unknown"}else{"Male"}})
  return(gender)
}

filter_affy <- function(x, mthr=4, mcv=NULL){
  pn <- rownames(x)
  gn <- affy2sym(pn)[,2]
  
  filt <- gn!='---' & !grepl('AFFX-', pn)
  
  if(!is.null(mthr)){
    filt <- filt & rowMeans(x) >= mthr
  }
  
  if(!is.null(mcv)){
    #require(genefilter)
    xcv <- rowSds(x)/rowMeans(x)
    filt <- filt & xcv >= mcv
  }
  
  return(which(filt))
}

ensembl2sym <- function(eid, concatenate=F, striptail=T, transcript=F){
  data(ensembl38)
  #ens <- readxl::read_excel('~/work/data/reference_genomes/hg38/ensembl_biomart_GRCh38_gene_info.xlsx')
  #ens <- readxl::read_excel('data/ensembl_biomart_GRCh38_gene_info.xlsx')
  
  ens <- ensembl38 %>% dplyr::select(`Gene stable ID`, `Gene name`) %>% unique
  if(transcript) ens <- ensembl38 %>% dplyr::select(`Transcript stable ID`, `Gene name`) %>% unique

  if(striptail) eid <- str_replace(eid, "\\.\\d+", "")
  ensid <- match(eid, ens$`Gene stable ID`)
  deprecated_ensid <- which(is.na(ensid))
  ercc <- grep("ERCC", eid)
  multimap_ensid <- which(ens$`Gene name`[ensid] %in% names(which(table( ens$`Gene name`[ensid])>1)))
  #unique_sym <- unique(ens$`Gene name`[ensid])
  #uid <- which(ens$`Gene name`[ensid] %in% unique_sym)
  
  if(concatenate) {
    sym <- paste0(ens$`Gene name`[ensid], "_", str_sub(eid, 10,15))
  } else{
    sym <- ens$`Gene name`[ensid]
  }
  return(list(genename=sym, deprecated_ensid=deprecated_ensid, ercc=ercc, multimap_ensid=multimap_ensid))
}

GSA.read.gmt=function(filename){
  #
  ## Read in and parse a gmt file (gene set file) from the  Broad institute
  # this is tricky, because each lines (geneset) has a variable length
  #  I read the file twice, first to pick up the geneset name and description
  # in the   first two  columns, then I read it all in as a long string
  
  # The beginning and end of each gene set in the string
  # is determined by matching
  # BOTH on  geneset name and description (since geneset names sometimes
  # occur as genenames elsewhere in the file)
  
  a=scan(filename,what=list("",""),sep="\t", quote=NULL, fill=T, flush=T,multi.line=F)
  geneset.names=a[1][[1]]
  
  geneset.descriptions=a[2][[1]]
  
  dd=scan(filename,what="",sep="\t", quote=NULL)
  
  
  nn=length(geneset.names)
  n=length(dd)
  ox=rep(NA,nn)
  
  ii=1
  for(i in 1:nn){
    cat(i)
    while((dd[ii]!=geneset.names[i]) | (dd[ii+1]!=geneset.descriptions[i]) ){
      ii=ii+1
    }
    ox[i]=ii
    ii=ii+1
  }
  
  genesets=vector("list",nn)
  
  for(i in 1:(nn-1)){
    cat(i,fill=T)
    i1=ox[i]+2
    i2=ox[i+1]-1
    geneset.descriptions[i]=dd[ox[i]+1]
    genesets[[i]]=dd[i1:i2]
  }
  
  geneset.descriptions[nn]=dd[ox[nn]+1]
  genesets[[nn]]=dd[(ox[nn]+2):n]
  out=list(genesets=genesets,geneset.names=geneset.names, geneset.descriptions=geneset.descriptions)
  class(out)="GSA.genesets"
  return(out)
}

library(Biobase)
setClass("ExpressionSetDC", contains="ExpressionSet", representation(fileinfo='list'))
setMethod("initialize", "ExpressionSetDC", 
          function(.Object, ...){
            callNextMethod(.Object, ...)
          })


getSurv <- function(clin, outcome='relapse', output='Surv'){
  require(survival)
  if(output=='Surv'){
    if(outcome=='relapse') {
      return(Surv(clin$pmwg.t_frelapse, clin$pmwg.ind_relapse=="Yes"))
    } else {
      return(Surv(clin[[paste0("t_", outcome)]], clin[[paste0(outcome, "_event")]]==1))
    }
  } else if(output=='df'){
    if(outcome=='relapse') {
      return(data.frame(time=clin$pmwg.t_frelapse, status=ifelse(clin$pmwg.ind_relapse=="Yes",1,0)))
    } else {
      return(data.frame(time=clin[[paste0("t_", outcome)]], status=ifelse(clin[[paste0(outcome, "_event")]]==1,1,0)))
    }
  }
}

mysurvfn <- function(vec, clinfile, param, thr=0.2, title=NULL){
  require(survminer)
  filt <- slicedat(vec, thr=thr)
  tempdf <- getSurv(clinfile[filt$id, ], param, output='df')
  tempdf$label <-filt$label
  ggsurvplot(survfit(Surv(time,status)~factor(label), data=tempdf), data=tempdf, pval=TRUE, conf.int = TRUE, title=title)[[1]]
}


plotSurvCov <- function(clin, outcome='relapse', covariate=NULL, ylim=c(0,1), main=NULL, 
                        print=T, returnmod=F, pval=F){
  s <- getSurv(clin, outcome)
  M=main
  if(is.null(M)) M=outcome
  mod <- getSurvReg(clin, outcome, covariate)
  if(pval) {pv <- round(mod$coefficients[5],4)} else {pv=''} 
  plot(survfit(s~covariate), col=1:length(unique(covariate)), 
       ylim=ylim, main=paste0(M, '\n',pv), ylab='Fraction', xlab='Time')
  
  if(print) print(mod)
  if(returnmod) return(mod)
}

getSurvReg <- function(clin, outcome='relapse', covariate=NULL, ylim=c(0,1)){
  s <- getSurv(clin, outcome)
  summary(coxph(s~factor(covariate)))
}


GSscore <- function(gset, dat, B=100){
  probes <- rownames(dat)
  gset2 <- update_gmt_obj(gset, probes)
  setscore <- list()
  
  for(f in 1:length(gset2$genesets)){
    gs <- gset2$geneset.names[f]
    print(gs)
    gsp <- gset2$genesets[[f]]
    np <- length(gsp)
    
    gsd <- colMeans(dat[gsp, ])
    gsd_r <- c()
    for(B in 1:100){
      gsd_r <- rbind(gsd_r, colMeans(dat[sample(probes, np),]))
    }
    tempscore <- c()
    for(g in 1:ncol(gsd_r)){
      tempscore <- c(tempscore, (gsd[g]-mean(gsd_r[,g]))/sd(gsd_r[,g]))
    }
    setscore[[f]] <- tempscore  
    
  }
  setscoremat <- ldply(setscore)
  rownames(setscoremat) <- gset2$geneset.names
  return(setscoremat)
}

surv <- function(m, tmax=NULL){
  if(is.null(tmax)) tmax <- max(m)
  out <- data.frame(t = 1:tmax, p = sapply(1:tmax, function(x) 100*sum(m >= x)/length(m)))
  return(out)
}

setcompare <- function(A, B, lengths=F){
  result <- list()
  result$uqA <- setdiff(A,B)
  result$uqB <- setdiff(B,A)
  result$common <- intersect(A,B)
  if(lengths){
    return(sapply(result, length))
  } else{
    return(result)
  }
}

allequal <- function(x){
  if(length(unique(x)) == 1){
    return(TRUE)
  } else {return(FALSE)}
}

madz <- function(x){
  z <- abs(x-median(x))/mad(x)
  return(z)
}

boxplot2 <- function(f, ...){
  pv <- anova(lm(f))$Pr[1]
  boxplot(f, main=round(pv, 5), pch=19, ...)
}

###################
# CompBio tools
###################
ProbeSelect <- function(cdat, sdat, zcutoff=1.5){
  cmeans <- matrix(rep(rowMeans(cdat), ncol(sdat)), ncol=ncol(sdat), byrow=F)
  csd <- matrix(rep(apply(cdat,1,sd), ncol(sdat)), ncol=ncol(sdat), byrow=F)
  zscores <- (sdat-cmeans)/csd
  
  N <- ncol(sdat)
  q <- 1-pnorm(zcutoff)
  mup <- apply(zscores, 1, function(x) sum(x>=zcutoff))
  pbinomup <- sapply(mup, function(x) pbinom(x, N, q))
  qbinomup <- p.adjust(pbinomup, 'fdr')
  
  mdn <- apply(zscores, 1, function(x) sum(x <= -zcutoff))
  pbinomdn <- sapply(mdn, function(x) pbinom(x, N, q))
  qbinomdn <- p.adjust(pbinomdn, 'fdr')
  
  # Combine q-values to get lesser value
  q <- c()
  for(i in 1:length(qbinomup)){
    q[i] <- min(qbinomup[i], qbinomdn[i])
  }
  return(q)
}

selectProbes <- function(exprs.df,ctrl_mu,ctrl_sig,pval_cut=1e-5,absZcut=1.5) {
  # Author: Raghavendra Hosur et al.
  # /home/rhosur/Projects/PersonalizedMed/2015/probeselect.R. Sourced on 09/30/2015
  #exprs.df is in the form of probes-by-patients
  #first makes sure everything is in the same order -- assuming ctrl_mu and ctrl_sig are in the same order
  comm_probes = names(ctrl_mu)[which(names(ctrl_mu)%in% rownames(exprs.df))] 
  print(paste("Length of common probes:",length(comm_probes)))
  exprs.df=exprs.df[comm_probes,]
  print(exprs.df[1:5,1:5])
  ctrl_mu = ctrl_mu[comm_probes]
  ctrl_sig = ctrl_sig[comm_probes]
  print(ctrl_mu[1:5]) 
  print(ctrl_sig[1:5])
  num_pats = dim(exprs.df)[2]
  exprs.df.zscores = apply(exprs.df,2,FUN=function(x){(x-ctrl_mu)/ctrl_sig})
  print(exprs.df.zscores[1:5,1:5])  
  exprs.df.succ = apply(exprs.df.zscores,1,FUN=function(x){sum(abs(x)>=absZcut)})
  p_select = pnorm(absZcut,lower.tail=F)+pnorm(-1.0*absZcut,lower.tail=T)
  exprs.df.binomP = lapply(exprs.df.succ,FUN=function(x){sum(dbinom(x:num_pats,num_pats,p_select))})
  probes_selected = names(which(exprs.df.binomP <= pval_cut))   
  return(probes_selected)
}

selProbes <- function(pdat, case, control, fdr=0.05, de=1.2, varratio=0.75, minexp=5, absZcut=1.5, plot=F, method='all'){
  prn <- rownames(pdat)
  # Probes with no annotations
  f00 <- which(affy2sym(rownames(pdat))[,2] != '---')
  prn00 <- prn[f00]
  
  # Minimum expr
  f0 <- which(rowMeans(pdat[f00, ]) >= minexp)
  prn0 <- prn00[f0]; f0 <- f00[f0]
  
  # DE
  #require(limma)
  design <- rep('CONTROL', ncol(pdat))
  design[case] <- 'CASE'
  fit <- eBayes(lmFit(pdat[f0,], model.matrix(~design)))
  prn1 <- rownames(topTable(fit, p.value=fdr, lfc=log2(de), num=Inf))
  if(tolower(method)=='de') return(prn1)
  
  # Var test
  casevar <- apply(pdat[f0,case],1,var)
  ctrvar <- apply(pdat[f0,control],1,var)
  rstat <- ctrvar/casevar
  fstat <- sapply(rstat, function(x) pf(x, length(control)-1, length(case)-1, lower.tail=T))
  prn2 <- names(which(p.adjust(fstat, 'fdr') < fdr & rstat <= varratio))
  if(tolower(method)=='var') return(prn2)
  
  # ProbeSelect
  ctrlm <- rowMeans(pdat[f0, control])
  ctrlsd <- apply(pdat[f0, control],1,sd)
  zscores <- t(sapply(1:length(ctrlm), function(i) (pdat[f0[i],case] - ctrlm[i])/ctrlsd[i]))
  exprs.df.succ <- apply(zscores, 1, function(x) sum(abs(x) >= 1.5))
  p_select = pnorm(absZcut,lower.tail=F)+pnorm(-1.0*absZcut,lower.tail=T)
  num_pats <- length(case)
  exprs.df.binomP = lapply(exprs.df.succ,FUN=function(x){sum(dbinom(x:num_pats,num_pats,p_select))})
  prn3 = prn0[which(p.adjust(exprs.df.binomP, 'fdr') < fdr)]
  if(tolower(method)=='ps') return(prn3)
  
  prn <- unique(c(prn1, prn2, prn3))
  if(plot) set_intersect <- venn(list(DEG=prn1, ProbeSelect=prn3, VarDiff=prn2))
  if(tolower(method)=='all') return(prn)
}

intergrep <- function(x, term, diff=FALSE){
  if(diff) return(x[-grep(term,x)])
  return(x[grep(term,x)])
}

melt.list2 <- function(x){
  xl <- lapply(x, melt)
  ox <- c()
  for(i in 1:length(xl)){
    xi <- xl[[i]]
    xi <- cbind(xi, names(xl)[i])
    xi$id <- rownames(xi)
    colnames(xi)[2] <- 'var'
    ox <- rbind(ox, xi)
  }
  rownames(ox) <- NULL
  return(ox)
}


limma2mat <- function(x,y,labels=NULL,fdr=0.05,thr=log2(1.5),mfilt=NULL,cfilt=TRUE){
  require(limma)
  if(!is.null(labels)){
    des <- model.matrix(~labels)
  } else {
    des <- model.matrix(~c(rep("GroupA", ncol(x)), rep("GroupB", ncol(y))))
  }
  
  # Check if rownames match
  if(!all(rownames(x)==rownames(y))) stop('Rownames don\'t match')
  
  dmat <- cbind(x,y)
  
  # Mean filter
  if(!is.null(mfilt)){
    dmat <- dmat[rowMeans(dmat)>=mfilt, ]
  }
  
  # Control probe filter (Affy only)
  if(cfilt & sum(grepl('AFFX', rownames(dmat)))>0) dmat <- dmat[-grep('AFFX-',rownames(dmat)), ]
  
  out <- list()
  out$fit <- eBayes(lmFit(dmat, des))
  out$ttable <- topTable(out$fit,p.value=  fdr, num = Inf,lfc = thr)
  out$res <- summary(decideTests(out$fit, lfc = thr, p.value = fdr))
  out$param <- data.frame(fdr=fdr, thr=thr, rowmean_filt=ifelse(is.null(mfilt),"NULL", mfilt), 
                          control_probe_filt=ifelse(is.null(cfilt),"NULL", cfilt))
  return(out)
}

limma1mat <- function(x,labels=NULL,fdr=0.05,thr=log2(1.5), mfilt=NULL, cfilt=TRUE){
  require(limma)
  if(is.null(labels)) stop('Must provide labels for limma 1 mat')
  des <- model.matrix(~labels)
  
  # NA filter
  naf <- which(is.na(labels))
  if(length(naf) > 0){
    print(paste("Removing NA samples ", naf, sep=" "))
    x <- x[,-naf]
    labels <- labels[-naf]
  }
  # Mean filter
  if(!is.null(mfilt)){
    x <- x[rowMeans(x)>=mfilt, ]
  }
  # Control probe filter (Affy only)
  if(cfilt & sum(grepl('AFFX', rownames(x)))>0) x <- x[-grep('AFFX-',rownames(x)), ]
  
  out <- list()
  out$fit <- eBayes(lmFit(x, des))
  out$ttable <- topTable(out$fit,p.value=  fdr, num = Inf,lfc = thr)
  out$res <- summary(decideTests(out$fit, lfc = thr, p.value = fdr))
  out$param <- data.frame(fdr=fdr, thr=thr, rowmean_filt=ifelse(is.null(mfilt),"NULL", mfilt), 
                          control_probe_filt=ifelse(is.null(cfilt),"NULL", cfilt))
  return(out)
}

limma1matfrac <- function(x,labels=NULL,fdr=0.05,thr=log2(1.5), mfilt=NULL, cfilt=TRUE){
  # Perform test on only a fraction of samples. 
  if(is.null(labels)) stop("Labels required for limma1matfrac")
  require(limma)
  sfilt <- which(!is.na(labels))
  x <- x[,sfilt]
  labels <- labels[sfilt]
  
  des <- model.matrix(~labels)
  
  # Mean filter
  if(!is.null(mfilt)){
    x <- x[rowMeans(x)>=mfilt, ]
  }
  # Control probe filter (Affy only)
  if(cfilt & sum(grepl('AFFX', rownames(x)))>0) x <- x[-grep('AFFX-',rownames(x)), ]
  
  out <- list()
  out$fit <- eBayes(lmFit(x, des))
  out$ttable <- topTable(out$fit,p.value=  fdr, num = Inf,lfc = thr)
  out$res <- summary(decideTests(out$fit, lfc = thr, p.value = fdr))
  out$param <- data.frame(fdr=fdr, thr=thr, rowmean_filt=ifelse(is.null(mfilt),"NULL", mfilt), 
                          control_probe_filt=ifelse(is.null(cfilt),"NULL", cfilt))
  return(out)
}

limmacov <- function(x,covariates=NULL,fdr=0.05,thr=log2(1.5), mfilt=NULL, cfilt=TRUE, impute=FALSE){
  require(limma)
  if(is.null(covariates)) stop('Must provide covariates data frame')
  
  if(impute) {
    cm <- apply(covariates, 2, function(z) ifelse(is.numeric(z),mean(z),"CHR"))
    for(i in 1:ncol(covariates)) covariates[which(is.na(covariates)),i] <- cm[i]
  }
  des <- model.matrix(~., covariates)
  
  # Mean filter
  if(!is.null(mfilt)){
    x <- x[rowMeans(x)>=mfilt, ]
  }
  # Control probe filter (Affy only)
  if(cfilt & sum(grepl('AFFX', rownames(x)))>0) x <- x[-grep('AFFX-',rownames(x)), ]
  
  out <- list()
  out$fit <- eBayes(lmFit(x, des))
  out$ttable <- topTable(out$fit,p.value=  fdr, num = Inf,lfc = thr)
  out$res <- summary(decideTests(out$fit, lfc = thr, p.value = fdr))
  out$param <- data.frame(fdr=fdr, thr=thr, rowmean_filt=ifelse(is.null(mfilt),"NULL", mfilt), 
                          control_probe_filt=ifelse(is.null(cfilt),"NULL", cfilt))
  return(out)
}

remove_duplicate_rows <- function(mat, cname="Gene", method=c('average', 'remove'), sep=NA, output='data.frame'){
  if(length(method) > 0) method='average'
  if(!"data.frame" %in% class(mat)) mat <- data.frame(mat, stringsAsFactors = F)
  cid <- match(cname, colnames(mat))
  if(!is.na(sep)) mat[[cname]] <- str_split(mat[[cname]], pattern = sep, simplify = T)[,1]
  dups <- names(which(table(mat[[cname]])>1))
  if(length(dups) == 0) {
    outmat <- mat[,-cname]
    rownames(outmat) <- mat[[cname]]
    return(outmat)
  }
  outmat <- mat[-match(dups, mat[[cname]]),-cid]
  
  for(i in 1:length(dups)){
    didx <- which(mat[[cname]]==dups[i])
    if(method=='remove'){
      vec <- mat[didx[1],-cid]
    } else{
      vec <- colMeans(mat[didx,-cid])
    }
    outmat <- rbind(outmat, vec)
    rownames(outmat)[nrow(outmat)] <- dups[i]
  }
  if(output=='matrix') outmat <- data.matrix(outmat)
  return(outmat)
}

format_results_table <- function(tt, type='deseq'){
  if(type=='deseq') {
    out <- tt %>% dplyr::select(1,2,3,7)
    colnames(out) <- c('Gene', 'AvgExpr', 'log2FoldChange', 'FDR')
  }
  if(type=='limma') {
    out <- tt %>% tbl_df %>% rownames_to_column('row') %>% dplyr::select(1,3,2,6)
    colnames(out) <- c('Gene', 'AvgExpr', 'log2FoldChange', 'FDR')
  }
  return(out)
}

topTable2 <- function(x, ...){
  t <- topTable(x, ...)
  t <- cbind(t, affy2sym(rownames(t)))
  return(t)
}

volcanoDESeq <- function(rmat, gstrip=F, top=10, vthr=1.5, fdr=0.05, print=TRUE, title=NULL){
  #if(gstrip) rmat$Gene <- str_split(rmat$row, "_", simplify=T)[,1]
  rmat <- data.frame(rmat)
  rmat$col <- ifelse(rmat$padj < fdr & abs(rmat$log2FoldChange) > log2(vthr), 'red', 'lightgray')
  out <- ggplot(rmat, aes(log2FoldChange, -log10(padj))) + geom_point(col=rmat$col) + theme_bw() + xlab('log2FoldChange') + 
    ylab('-log10(FDR)') + geom_vline(xintercept=c(log2(vthr), log2(1/vthr))) + geom_hline(yintercept=-log10(fdr))
  if(!is.null(title)) out <- out + ggtitle(title)
  if(print) {
    print(out)
  } else{
    return(out)
  }
}


fitEnsemble <- function(trdat, cl, method='all'){
  outli <- list()
  # if(method=='all'){
  #   # PAM
  #   #######################
  #   print('Fitting PAMR')
  #   require(pamr)
  #   outli$pamr.train <- pamr.train(list(x=trdat, y=factor(cl)), n.threshold=100)
  #   outli$pamcv <- pamr.cv(outli$pamr.train, list(x=trdat, y=factor(cl)))
  # }
  if(method=='all'){
    print('Fitting Lasso')
    require(glmnet)
    # Lasso
    #######################
    FAM <- ifelse(length(unique(cl))==2, 'binomial', 'multinomial')
    outli$lasso.cvfit <- cv.glmnet(x=t(trdat), y=factor(cl), family=FAM, standardize=T, alpha=0.5)
  }
  if(method=='all'){
    print('Fitting RandomForest')
    require(randomForest)
    # Random Forest
    #######################
    outli$rf.fit <- randomForest(x=t(trdat), y=factor(cl), importance = T)
  }
  if(method=='all'){
    print('Fitting SVM')
    require(e1071)
    # SVM
    #######################
    outli$svm.fit <- tune.svm(x=t(trdat), y=factor(cl))
  }
  if(method=='all'){
    print('Fitting C5.0')
    require(C50)
    # C50
    #######################
    outli$C50.fit <- C5.0(x=t(trdat), y=factor(cl))
  }
  if(method=='all'){
    print('Fitting XGBOOST')
    require(xgboost)
    xtr <- xgb.DMatrix(t(trdat), label=cl)
    FAM <- ifelse(length(unique(cl))==2, "binary:logistic", "multi:softmax")
    NUMCLASS = length(unique(cl))
    PLIST <- list(objective = FAM, max.depth =3, eta = 1, nthread=2)
    if(NUMCLASS > 2) PLIST$num_class=NUMCLASS
    outli$xgb.fit <- xgboost(data=xtr, nrounds=3, nfold=5, params=PLIST)
  }
  
  # if(method='all'){
  #   require(gbm)
  #   print('Fitting GBM')
  #   # GBM
  #   #######################
  #   temp <- list(x=t(trdat), y=factor(cl))
  #   outli$gbm.fit <- gbm(y~x, temp, distribution='multinomial', train.fraction=0.8, n.trees=1000, shrinkage = 0.05, cv.folds=10)
  # }
  # if(method='all'){
  #   require(e1071)
  #   print('Fitting NaiveBayes')
  #   # NaiveBayes
  #   #######################
  #   temp <- data.frame(cbind(t(trdat), cl))
  #   outli$nb.fit <- gbm(cl~., temp)
  # }
  # if(method='all'){
  #   require(klaR)
  #   print('Fitting RDA')
  #   # RDA
  #   #######################
  #   temp <- list(x=t(trdat), y=factor(cl))
  #   outli$rda.fit <- rda(y~x, temp,  gamma = 0.05, lambda = 0.2))
  # }
  # if(method='all'){
  #   require(pls)
  #   print('Fitting PLS')
  #   # PLS
  #   #######################
  #   temp <- list(x=t(trdat), y=cl)
  #   outli$mvr.fit <- mvr(y~x, temp, validation="CV")
  # }
  return(outli)
}


predEnsembl <- function(trm, ndat, plot=F, aggregate=c('average', 'weighted')[1]){
  outli <- list()
  runmodels <- names(trm)
  
  # # PAMR
  # if('pamcv' %in% runmodels){
  #   thr <- max(0, trm$pamcv$threshold[which.min(trm$pamcv$error)], na.rm=T)
  #   outli$pamr <- as.numeric(pamr.predict(trm$pamr.train, newx=ndat, type="class", threshold=thr))
  # }
  
  # Lasso
  if('lasso.cvfit' %in% runmodels){
    outli$lasso <- as.numeric(predict.cv.glmnet(trm$lasso.cvfit, newx=t(ndat), 
                                                type='class', s='lambda.min'))
  }
  
  # RandomForest
  if('rf.fit' %in% runmodels){
    outli$randomforest <- as.numeric(predict(trm$rf.fit, newdata=t(ndat))) - 1
  }
  
  # SVM
  if('svm.fit' %in% runmodels){
    outli$svm <- as.numeric(predict(trm$svm.fit$best.model, newdata=t(ndat))) - 1
  }
  
  # C5.0
  if('C50.fit' %in% runmodels){
    outli$C5.0 <- as.numeric(predict(trm$C50.fit, newdata=t(ndat))) - 1
  }
  
  # XGBOOST
  if('xgb.fit' %in% runmodels){
    xp <- as.numeric(predict(trm$xgb.fit, newdata=t(ndat)))
    if(length(unique(xp)) > 5) xp <- ifelse(xp > 0.5, 1, 0)
    outli$xgboost <- xp
  }
  
  outdf <- list()
  outdf$preddf <- as.data.frame(outli)
  
  if(aggregate=='average'){
    outdf$predmax <- as.numeric(apply(outdf$preddf, 1, function(x) names(sort(-table(x)))[1]))
  } else if(aggregate=='weighted'){
    w <- apply(outdf$preddf, 2, function(x) dist(t(data.frame(x=x, y=cl)), method='manhattan'))
    w <- 1-w/length(cl)
    outdf$predmax <- round(as.vector(data.matrix(outdf$preddf) %*% w)/ncol(outdf$preddf))
  }
  names(outdf$predmax) <- colnames(ndat)
  
  if(plot) {
    plot(hclust(dist(t(outdf$preddf), 'manhattan'), method='ward.D2'))
    require(gplots)
    for(cl in 1:max(as.numeric(apply(outdf$preddf,2,max)))){
      par(ask=TRUE) 
      venn(apply(outdf$preddf, 2, function(x) which(x==cl)))
    }
  }
  
  return(outdf)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

figcapfn = local({
  i = 0
  function(x) {
    i <<- i + 1
    paste('Figure ', i, ': ', x, sep = '')
  }
})

fill_missing <- function(x){
  x <- as.numeric(x)
  id <- which(is.na(x))
  x[id] <- mean(x[-id])
  return(x)
}

biogenpalette <- c('#1c5a7d', '#578196', '#2573ba', '#6dad46', '#679acb',
                   '#519643', '#7cc3e2', '#99ca3c', '#c7dd72', '#dde5ae', '#7c878e')


getGlmnetCoef <- function(mod, s='lambda.min', intercept=F){
  c <- coef.cv.glmnet(mod, s=s)
  o <- abs(c[which(as.vector(c)!=0),]) %>% sort(decreasing = T) %>% as.data.frame
  o$var <- rownames(o)
  colnames(o)[1] <- 'coef'
  if(!intercept) o <- o[-match("(Intercept)", o$var),]
  return(o)
}

getGlmnetTopN <- function(mod, N=10, intercept=F){
  S <- mod$lambda[which(mod$nzero > N)[1]]
  return(getGlmnetCoef(mod, s=S, intercept = intercept))
}

getTrnProbes <- function(dat, cov){
  pid <- apply(dat, 1, function(x) summary(lm(x ~ cov))$coef[2, 4])   %>% 
    Filter(function(x) x<=0.1,.)   %>% names %>% match(., rownames(dat))
  rpid <- resid(lm(t(dat[pid,])~cov)) %>% t
  rpidcor <- cor(t(rpid)) %>% abs
  rpidcorlist <- data.frame(con =  rpidcor %>% apply(1, function(x) length(which(x>.5))), var = apply(rpid, 1, var))
  rpidcorlist <- rpidcorlist[order(rpidcorlist$con, rpidcorlist$var, decreasing = T), ]
  rmlist <- c()
  if(rpidcorlist$con %>% max > 1){
    for(i in rownames(subset(rpidcorlist, con > 1))){
      if(! i %in% rmlist) rmlist <- c(rmlist, setdiff(names(which(rpidcor[i, ] >= 0.5)), i))
    }
  }
  return(setdiff(rownames(dat)[pid], rmlist))
}

plotCO <- function(clin, cl, col=NULL){
  par(mfrow=c(4,4))
  for(F in paste0('pmwg.', c('GDVOLBL', 'GDLESBL', 'T2VOLBL',  'T1VOLBL', 'PGLB_MSSS', 'NBVBL', 
             'newT2_24wk', 'volGD_24wk', 'cntGD_24wk', 'WBMTRBL', 'BVpch_24wk', 'ONSYRS', 
             'num_relapse', 'age', 'MSFCBL', 'BASEALCimp'))){
    filt <- rep(TRUE, nrow(clin))
    ymax = quantile(clin[[F]][filt], 0.9, na.rm=T)
    ymin = quantile(clin[[F]][filt], 0.05, na.rm=T)
    boxplot(clin[[F]][filt]~cl[filt], col=col,
            main=paste0(F, ': ', round(kruskal.test(clin[[F]][filt]~cl[filt])$p.value,3)), 
            ylim=c(ymin,ymax))
  }
}

plotCO2 <- function(clin, cl){
  newdf <- list()
  newdf$frac_prog6mo[['0']] <- mean(clin$pwmg.t_prog12w < 183)
  newdf$frac_edss6mo[['0']] <- mean(clin$pmwg.t_edssplus < 183)
  newdf$frac_relapse6mo[['0']] <- mean(clin$pmwg.t_frelapse < 183)
  for(K in 1:max(cl)){
    newdf$frac_prog6mo[[as.character(K)]] <- mean(clin$pwmg.t_prog12w[cl==K] < 183)
    newdf$frac_edss6mo[[as.character(K)]] <- mean(clin$pwmg.t_edssplus[cl==K] < 183)
    newdf$frac_relapse6mo[[as.character(K)]] <- mean(clin$pwmg.t_frelapse[cl==K] < 183)
  }
  return(newdf)
}



###################
# Cluster tools
###################

swapcluster <- function(cl, l, r){
  #L <==> R
  t <- cl
  t[cl==l] <- r
  t[cl==r] <- l
  return(t)
}

colLab <- function(n, labelColors, clusMember) {
  if(is.leaf(n)) {
    a <- attributes(n)
    # clusMember - a vector designating leaf grouping
    # labelColors - a vector of colors for the above grouping
    labCol <- labelColors[clusMember[match(a$label,names(clusMember))]]
    attr(n, "edgePar") <- c(a$nodePar, list(col = labCol, lab.cex=0.1))
  }
  n
}

paintdgram <- function(hc, clusmem, col=NULL, title=NULL){
  require(RColorBrewer)
  if(is.null(col)){
    labelColors <- brewer.pal(8, 'Set1')[1:max(clusmem)]
  } else{
    labelColors <- col
  }
  dhc <- as.dendrogram(hc)
  clusMember=clusmem
  dL <- dendrapply(dhc, colLab, labelColors, clusMember)
  op <- par(mar = c(2,4,4,2) + 0.1)
  t='Cluster'
  if(!is.null(title)) t = title
  plot(dL,  xlab="", sub="", xaxt="n", main=t)
  par(op)
}

pacc <- function(h, c, B=10000){
  coph <- cophenetic(h)
  d <- as.matrix(coph)
  #c <- c[match(attr(coph, 'Labels'), names(c))]
  cm <- match(c, attr(coph, 'Labels'))
  do <- c()
  db <- c()
  for(boot in 1:B){
    os <- sample(cm, replace = T)
    do[boot] <- mean(d[os,os])
    bs <- sample(nrow(d), length(cm))
    db[boot] <- mean(d[bs, bs])
  }
  z <- (mean(do)-mean(db))/(sd(db)/sqrt(B))
  return(list(zscore=z, pvalue=pnorm(z)))
}

match0 <- function(x,y){
  out <- match(x,y)
  out <- out[!is.na(out)]
  return(out)
}

median0 <- function(x) median(x, na.rm=T)
sum0 <- function(x) sum(x, na.rm=T)
mean0 <- function(x) mean(x, na.rm=T)

cut_number2 <- function(x, Q=4, dec=F, prefix="Q", value=NULL){
  if(!is.null(value)) {
    y <- ifelse(x >= value, "hi", "lo")
    return(y)
  }
  
  y <- ggplot2::cut_number(x, Q)
  l <- levels(y)
  if(dec) l <- rev(l)
  q <- paste0(prefix, 1:Q)
  return(q[match(y,l)])
}
