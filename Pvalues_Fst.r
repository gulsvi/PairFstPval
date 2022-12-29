
############################################
###       pairwise Fst and p-values      ###
############################################

Pval_Fst <- function(x,path,perm=10000){
## x : class ‘genind’  adegenet R-package
fst <- pairwise.fst(x, res.type="matrix")          
mat.perm <- lapply(1:perm, function(i) pairwise.fst(x, pop=sample(pop(x)), res.type="matrix"))
mean(c(fst[1,2] < na.omit(sapply(1:perm, function(i) mat.perm[[i]][1,2])), TRUE))
####
test12 <- as.randtest(na.omit(sapply(1:perm, function(i) mat.perm[[i]][1,2])), fst[1,2], alter="greater")
####
z <- list()
tab_pval<-matrix(NaN,dim(fst)[1],dim(fst)[2])
rownames(tab_pval)<-rownames(fst); colnames(tab_pval)<-colnames(fst)

 for(i in 1:(nrow(fst))){
   for(j in 1:nrow(fst)){
   sss<-as.randtest(na.omit(sapply(1:perm, function(k) mat.perm[[k]][i,j])), fst[i,j], alter="greater")
   z[[paste(rownames(fst)[i],rownames(fst)[j],sep="-")]] <-sss$pvalue
   if(i<j) {tab_pval[i,j]<-sss$pvalue} else {tab_pval[i,j]<-fst[i,j]}
   }
}
write.table(tab_pval,paste0(path,'pairwise_Fst_and_pvalues.txt'),sep='\t')
return(tab_pval)
}

