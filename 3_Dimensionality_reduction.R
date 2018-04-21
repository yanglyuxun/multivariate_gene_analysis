# update selected features from Python result
#X = as.matrix(read.csv('selectedX.csv',row.names=1,header=T))
#save(X,file='selectedX.rdata')
# load data
load('selectedX.rdata')

# Visualization --------------
library(ggcorrplot) # plot for corr
## corr matrix between patiens
corr = round(cor(t(X)), 1) 
ggcorrplot(unname(corr)) #;ggsave('corr_patients.png')
ggcorrplot(unname(ifelse(abs(corr)>0.5,corr,0)))
## corr matrix between microarrays
corr = round(cor(X), 1) 
ggcorrplot(unname(corr[1:100,1:100]))
corr0 = ifelse(abs(corr)>0.5,corr,0)
ggcorrplot(unname(corr0[1:100,1:100]))
ggcorrplot(unname(corr0[100:200,100:200]))
rm(corr,corr0)

# Dimensionality reduction -----------
## Dim Red method selection
library(dimRed)
libs = c('RANN', 'RSpectra', 'Rtsne', 'coRanking', 'diffusionMap',
         'energy','fastICA','ggplot2','graphics', 'igraph', 'kernlab', 'lle','loe',
         'optimx', 'pcaPP', 'rgl', 'scales', 'scatterplot3d', 'stats', 'testthat',
         'tidyr', 'vegan') # the libraries that "dimRed" need to use
for(l in libs){
  if(!require(l,character.only=T)){install.packages(l)}
}
embed_methods = dimRedMethodList() # all methods
quality_methods = dimRedQualityList() # all Quality measures
quality_results = matrix(NA, length(embed_methods), length(quality_methods),
                         dimnames = list(embed_methods, quality_methods))
embedded_data = list()
for (e in embed_methods) {
  message("embedding: ", e)
  embedded_data[[e]] = embed(X, e, .mute = character(0))
  for (q in quality_methods) {
    message("  quality: ", q)
    quality_results[e, q] = tryCatch(
      quality(embedded_data[[e]], q),
      error = function(e) NA
    )
    write.csv(quality_results,'tmp2.csv')
  }
}
save(embedded_data,file='embedded_data.rdata',compress=T,compression_level=9)
embed0 = embed(X, 'PCA')
embed0@data@data = as.matrix(read.csv('FAresult.csv'))
quality_results=rbind(quality_results,NA)
rownames(quality_results)[16]='FA'
for (q in quality_methods){
  quality_results['FA',q]=tryCatch(
    quality(embed0,q),error = function(e) NA);message(q)
}
write.csv(quality_results,'quality_results.csv')
save(quality_results,file='quality_results.rdata')
## rank the results (NOT selected)
# quality_order = quality_results[,1:7]+NA # initialization
# for(i in 1:ncol(quality_order)){
#   quality_order[,i]=rank(quality_results[,i])
# }
# sort(rowMeans(quality_order)) # the rank sort
## normalize the results
quality_normal = quality_results[,1:7]+NA # initialization
for(i in 1:ncol(quality_normal)){
  quality_normal[,i]=(quality_results[,i]-min(quality_results[,i]))/
    (max(quality_results[,i])-min(quality_results[,i]))
}
sort(rowMeans(quality_normal[,1:7]))
## plot the quality results
library(ggplot2)
plot_mat=function(m){
  df = data.frame(Index=NA,Value=NA,method=NA)
  p=1
  for(i in 1:nrow(m)){
    for(j in 1:(ncol(m)-1)){
      if(is.na(m[i,j])){next}
      df[p,1]=colnames(m)[j]
      df[p,2]=m[i,j]
      df[p,3]=rownames(m)[i]
      p=p+1
    }
  }
  df$method=factor(df$method)
  ggplot(df)+geom_line(aes(x=Index,y=Value,col=method,group=method))+
    geom_text(aes(label=method, col=method, x=Index,y=Value),alpha=0.7)+
    xlab('Quelity Index')+ylab('Quality Value')+theme_bw()+
  theme(axis.text.x = element_text(angle=20,hjust=1))
}
plot_mat(quality_results)

## plot the 2D graph for each method
GMresult = factor(as.matrix(read.csv('GaussianMixtureResult.csv')))
embedded_data[['FA']] = embedded_data[['PCA']] #initialize and read FA results
embedded_data[['FA']]@data@data = as.matrix(read.csv('FAresult.csv'))
plot_2d=function(data,main=NULL){ # a qplot wrapper
  #ggplot(NULL,aes(x=data[,1],y=data[,2],color=GMresult))+
  ggplot(NULL,aes(x=data[,1],y=data[,2]))+
    geom_point(alpha=0.5)+
    xlab(colnames(data)[1])+ylab(colnames(data)[2])+
    ggtitle(main)+
    theme_bw()
}
plot_result=function(name){ # a result plotting wrapper
  plot_2d(embedded_data[[name]]@data@data, main=name)
}
for(name in names(sort(rowMeans(quality_normal[,1:7])))){
  plot(plot_result(name))
}

## PCA analysis
library(pcaMethods) # use the package directly, rather than the wrapper
pc = pca(X,'svd',nPcs=nrow(X),completeObs=F)
#loadings = loadings(pc)
#scores = scores(pc)
slplot(pc, scoresLoadings=c(T, F))
sum(pc@R2[1:2])/sum(pc@R2) # 0.19969
plot_2d(pc@scores) # same as before
ggplot(NULL,aes(x=1:50,y=pc@R2[1:50]))+xlab('number of PC')+
  ylab('Variance')+geom_line()+geom_point()+theme_bw()
ggplot(NULL,aes(x=1:length(pc@R2cum),y=pc@R2cum))+geom_line()+geom_hline(yintercept=0.9)+theme_bw()
which_cum_p=function(p){
  #find the number of PCs that explain p of the total variance
  which(pc@R2cum>=p)[1]
}
which_cum_p(0.3) #4
which_cum_p(0.5) #13
which_cum_p(0.7) #34
which_cum_p(0.8) #54
which_cum_p(0.9) #91
which_cum_p(0.95) #124
which_cum_p(0.99) #184

## Kernel PCA
library(kernlab) # use the package directly, rather than the wrapper
kpc = kpca(X, kernel = "rbfdot", kpar = list(sigma = 1/ncol(X)))
plot(rotated(kpc),main=sum(kpc@eig[1:2])/sum(kpc@eig))
plot(pc@R2[1:50]/sum(pc@R2))
lines(pc@R2[1:50]/sum(pc@R2))
points(kpc@eig[1:50]/sum(kpc@eig))
lines(kpc@eig[1:50]/sum(kpc@eig))


## factor analysis 
# -->> Use Python (EM algorithm) because the following codes are too slow <<--
# library("psych")
# fafit = principal(X, nfactors=2, rotate="varimax")
# fafit = principal(X, nfactors=2, rotate="oblimin")
# fafit = fa(X, nfactors=2, rotate="oblimin",fm="wls") 
# fafit = fa(X, nfactors=2, rotate="varimax",fm="wls") 
# fafit = factanal(X,2)
# fafit = factanal(X,2,rotation="varimax")
FAX = as.matrix(read.csv('FAresult.csv')) #<-- read data created by python
plot(FAX)

