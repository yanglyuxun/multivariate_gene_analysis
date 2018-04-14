# Libraries -------------------
## install
# source("https://bioconductor.org/biocLite.R")
# biocLite("impute") 
# biocLite("pcaMethods") 
## import
library(impute) # for knn imputation
library(pcaMethods) # for other imputation methods

# Read data and reorganize -------------------
expr=read.csv(file="microarray_data_NEJM.csv",header=TRUE,stringsAsFactors=F)
LYMexpr=expr[,1:276]
response=read.csv(file="Patient_data_NEJM.csv",header=TRUE)
newresponse=response[order(response[,1]),]
idvec=NULL
for (j in 3:(dim(LYMexpr)[2]))
{
  idvec[j]=substr(unlist(strsplit(names(LYMexpr)[j],"[_]"))[2],4,6)
}
idvec[1:2]=names(LYMexpr)[1:2]
colnames(LYMexpr)=idvec
indicator=c(TRUE,TRUE,(as.numeric(idvec[3:276]))%in%response[,1])
LYMexpr4use=LYMexpr[,indicator]
newLYMexpr=LYMexpr4use[,c(1,2,order(names(LYMexpr4use)[3:242])+2)]
rm(expr,response,idvec,LYMexpr,indicator,LYMexpr4use,j)

# basic analysis -----------------
max(newLYMexpr[,3:ncol(newLYMexpr)],na.rm = TRUE)
summary(newLYMexpr[,3:ncol(newLYMexpr)])
summary(as.array(unlist(newLYMexpr[,3:ncol(newLYMexpr)])))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  -8.20   -0.39    0.00    0.01    0.40    8.37  182919 
plot(density(as.array(unlist(newLYMexpr[,3:ncol(newLYMexpr)])),na.rm=T))
plot(density(as.array(unlist(newLYMexpr[,6])),na.rm=T))

# Imputation Experiments ---------------------
## remove microarrays with too many missing values
X = t(newLYMexpr[,3:ncol(newLYMexpr)]) # the n*p data matrix
dim(X) # 240 7399
colnames(X) = unname(unlist(newLYMexpr['NAME']))
na_rate = function(x){
  ax = as.array(x)
  sum(is.na(x))/length(x)
}
total_na_rate = na_rate(X) # 10.30%
gene_na = rep(0,ncol(X))
for(i in 1:ncol(X)){gene_na[i] = na_rate(X[,i])}
hist(gene_na, xlab = 'NA rate', main='Histogram of NA rates')
sum(gene_na>=0.4)/length(gene_na)
# 4.85% of microarrays have >=40% missing values, !!!REMOVE in X!!!
X = X[,gene_na<0.4] 
dim(X) # 240 7040

## experiments on different imputation methods
impu_exp = function(X, r_na, n=100){
  # X= data matrix
  # r_na= the ratio of number of na to make
  # n= number of repeated sampling
  candidates = which(!is.na(X))
  NRMSE = function(X,impu_X,new_na){
    guess = impu_X[new_na]
    answer = X[new_na]
    sqrt(mean((guess-answer)^2)/var(answer))
  }
  methods = c('KNN2','KNN5','KNN10','LLS','PPCA','BPCA','SVD')
  ans = matrix(0,n,length(methods)) # the result matrix
  colnames(ans) = methods
  for(i in 1:n){
    new_na = sample(candidates, r_na*length(candidates))
    newX = X
    newX[new_na] = NA
    ans[i,'KNN2'] = NRMSE(X, t(impute.knn(t(newX),k=2)$data), new_na)
    cat(i,'KNN2',ans[i,'KNN2'],'\n');write.csv(ans,'tmp.csv')
    ans[i,'KNN5'] = NRMSE(X, t(impute.knn(t(newX),k=5)$data), new_na)
    cat(i,'KNN5',ans[i,'KNN5'],'\n');write.csv(ans,'tmp.csv')
    ans[i,'KNN10'] = NRMSE(X, t(impute.knn(t(newX),k=10)$data), new_na)
    cat(i,'KNN10',ans[i,'KNN10'],'\n');write.csv(ans,'tmp.csv')
    ans[i,'LLS'] = NRMSE(X, completeObs(llsImpute(newX,k=10,verbose=T,allVariables=T,maxSteps=10)), new_na)
    cat(i,'LLS',ans[i,'LLS'],'\n');write.csv(ans,'tmp.csv')
    ans[i,'PPCA'] = NRMSE(X, completeObs(pca(newX,"ppca")), new_na)
    cat(i,'PPCA',ans[i,'PPCA'],'\n');write.csv(ans,'tmp.csv')
    ans[i,'BPCA'] = NRMSE(X, completeObs(pca(newX,"bpca")), new_na)
    cat(i,'BPCA',ans[i,'BPCA'],'\n');write.csv(ans,'tmp.csv')
    ans[i,'SVD'] = NRMSE(X, completeObs(pca(newX,"svdImpute")), new_na)
    cat(i,'SVD', ans[i,'SVD'] ,'\n');write.csv(ans,'tmp.csv')
  }
  return(ans)
}
impu_comp = impu_exp(X,na_rate(X))
colMeans(impu_comp)
## plot the boxes
library(ggplot2)
box_mat=function(m){
  df  =data.frame(method=colnames(m),NRMSE=as.vector(t(m)))
  ggplot(df)+geom_boxplot(aes(x=method,y=NRMSE))+theme_bw()
}
box_mat(impu_comp)

# Imputation -----------
newX = completeObs(llsImpute(X,k=10,verbose=T,allVariables=T,maxSteps=10))
save(newX,file='imputedX.rdata')
write.csv(newX,'imputedX.csv')
