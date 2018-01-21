#################################################################
###data preprocessing
##log transformation for nonzero values for each feature column
##match case control study
##data include three parts-response-1,2
##paired group id
##features matrix 
##select 406 out of 472 features with missing percentage less than or equal to 0.05
##impute data by random samples from 5% tail region of assuming log normal distribution
#########################################################
setwd("/Users/yubingyao/Documents/Course materials/CS590V/Final project/data")
#response dataset
dat.resp<-read.table("ResponseVector.txt")
dim(dat.resp)
#feature dataset
feature.dat<-read.table("Data.txt")

feature.zero.id=apply(feature.dat,2,function(x){which(x==0)})
missing.perc=matrix(unlist(lapply(feature.zero.id,function(x){length(x)/dim(feature.dat)[1]})),c(1,dim(feature.dat)[2]))
length(which(missing.perc==0))
length(which(missing.perc<=0.05 & missing.perc>0))
length(which(missing.perc>0.05))

feature.nonzero<-feature.dat[,which(missing.perc==0)]
#select the features with mising percentage less than or equal to 5 percentage
feature.zerosel<-feature.dat[,which(missing.perc<=0.05 & missing.perc>0)]
#impute.tail.lognormal:assumiing each column of features follow log normal distributin
#with very small missing percentage, approximate assume mean of log normal distribution equal 
#to sample mean of log nonzero values, distribution standard deviation eqwual
#to sample standard deviation of log transformed nonzero values
#impute the zero values with random samples with tail region of log normal distribution
impute.tail.lognormal<-function(x){
 # x=feature.zerosel[,1]
  n.samp=sum(x==0)
  tail.prob=sum(x==0)/length(x)
  logmean=mean(log(x[x!=0]))
  logstd=sd(log(x[x!=0]))
  #tail.prob<-ifelse(tail.prob>0.5,1-tail.prob,tail.prob)
  val.thre<-min(c(qlnorm(tail.prob,meanlog=logmean,sdlog=logstd),x[x!=0]))
  samp.vec=rep(0,n.samp)
  for (i in 1:n.samp){
    sample.flag<-T
     while(sample.flag){
        samp<-rlnorm(1, meanlog = logmean, sdlog = logstd)
        if(samp<val.thre){
          samp.vec[i]=samp
           sample.flag=F
         }
     }
  }
  x[x==0]=samp.vec
  x
}
feature.zerosel.impute<-apply(feature.zerosel,2,impute.tail.lognormal)
feature.zerosel.impute.scale<-apply(feature.zerosel.impute,2,function(x){(log(x)-mean(log(x)))/sd(log(x))})
colnames(feature.zerosel.impute.scale)
#log transform the features, scale by deducting the mean divided by standard deviation
log.feature.nonzero=log(feature.nonzero)
feature.nonzero.scale<-apply(log.feature.nonzero,2,function(x){(x-mean(x))/sd(x)})
#final features after scale and impute and log transformation, selection of all the features
feature.transform<-cbind(feature.zerosel.impute.scale,feature.nonzero.scale)
#pair group id
dat.pairid<-read.table("MatchedPairIndicator.txt")
dim(dat.pairid)
unique(sort(as.matrix(dat.pairid)))
table(dat.pairid)
dat_tol=cbind(feature.dat,dat.resp,dat.pairid)
dim(dat_tol)
colnames(dat_tol)<-c(colnames(feature.dat),"response","pair_id")
##pca
##output with first 10 principal components
feature.pca <- prcomp(feature.transform,
                 center = TRUE,
                 scale. = TRUE) 
dim(feature.pca)
feature.transform[1,]
dim(feature.pca$rotation)
pca.feature<-feature.transform%*%feature.pca$rotation
? prcomp
##pls
##output with first 15,20 components
library(plsdepot)
library(pls)
pls_feature = plsreg1(feature.transform, dat.resp, comps = 15)
#15 components account for 95% correlation between the features and response
#20 components account for 99% correlation between the features and response
#10 components account for 88% correlation between the features and response
comp.15.pls<-pls_feature$x.scores

dim(comp.15.pls)
pls_feature2 = pls(feature.transform, dat.resp, comps = 15)
sum(pls_feature$R2)
dim(pls_20components)
#the pls transformed 15 components 
pls_15components<-pls_feature$x.scores
pls_feature_20 = plsreg1(feature.transform, dat.resp, comps = 20)
pls_20components<-pls_feature_20$x.scores

colnames(pls_20components)<-sapply(1:ncol(pls_20components),function(x){gsub('t','PLS_Component',colnames(pls_20components)[x])})
pls_20components_y<-cbind(pls_20components,diseased_group)
write.csv(pls_20components_y, file = "PLS_20components_diseasegroup.csv",row.names=FALSE)
##paired t test
#transform the features data grouped by pairid and response status(1,2)
#each feature column to 2 columns, first column-1-case diseased group, second 
#column healthy control group-2, with each row the observations with same
#paired id from 1 to 86

#paired t test
#select 10 or 20 features with smallest p values
feature.disease<-feature.dat[dat.resp==1,]
feature.health<-feature.dat[dat.resp==1,]
dat_tol[,1][dat.resp==2 & dat.pairid==1]
pair_dat_tol=sapply(unique(dat.pairid),function(x){apply(dat_tol,2,function(y){cbind(y[dat.resp==1 & dat.pairid==x],y[dat.resp==2 & dat.pairid==x])})})
mat_trans_pairid=matrix(NA,nrow=dim(unique(dat.pairid))[1],ncol=2*dim(feature.transform)[2])
for (i in unique(sort(as.matrix(dat.pairid)))){
  for (j in 1:dim(feature.transform)[2]){
    mat_trans_pairid[i,2*j]<-feature.transform[,j][dat.resp==2 & dat.pairid==i]
    mat_trans_pairid[i,2*j-1]<- feature.transform[,j][dat.resp==1 & dat.pairid==i]
  }
  }
  
pairtt_pvalue=sapply(1:dim(feature.transform)[2],function(i){t.test(x=mat_trans_pairid[,2*i-1],y=mat_trans_pairid[,2*i],paired=T)$p.value})
pairtt_pvalue[order(pairtt_pvalue)[1:20]]
feature.10lowestpvalue<-feature.transform[,order(pairtt_pvalue)[1:10]]
feature.20lowestpvalue<-feature.transform[,order(pairtt_pvalue)[1:20]]
dat.resp.new<-as.matrix(dat.resp.new)
diseased_group<-apply(dat.resp.new,1,function(x){ifelse(x==1,'Diseased','Healthy')})
diseased_group<-as.matrix(diseased_group)
colnames(diseased_group)<-"Disease_group"
feature.10lowestpvalue.diseasegroup<-cbind(feature.10lowestpvalue,diseased_group)
feature.20lowestpvalue.diseasegroup<-cbind(feature.20lowestpvalue,diseased_group)
colnames(feature.10lowestpvalue.diseasegroup)<-c("Metabolite322" , "Metabolite28" ,          "Metabolite117"  ,        "Metabolite132",
                            "Metabolite70" , "Metabolite169"  ,        "Metabolite215" ,         "Metabolite92" ,         
                            "Metabolite203",          "Metabolite29","Disease_group")
# Write CSV in R
colnames(feature.20lowestpvalue.diseasegroup)
colnames(feature.20lowestpvalue.diseasegroup)<-c(sapply(1:(length(colnames(feature.20lowestpvalue.diseasegroup))-1),function(x){paste("Metabolite",unlist(strsplit(colnames(feature.20lowestpvalue.diseasegroup)[x],'V'))[2],sep='')
}),"Disease_group")
write.csv(feature.10lowestpvalue.diseasegroup, file = "feature_10lowestpvalue_diseasegroup.csv",row.names=FALSE)
write.csv(feature.20lowestpvalue.diseasegroup, file = "feature_20lowestpvalue_diseasegroup.csv",row.names=FALSE)

colnames(feature.transform)<-sapply(1:(length(colnames(feature.transform))),function(x){paste("Metabolite",unlist(strsplit(colnames(feature.transform)[x],'V'))[2],sep='')
})
feature.transform.diseasegroup<-cbind(feature.transform,diseased_group)
write.csv(feature.transform.diseasegroup, file = "feature_diseasegroup.csv",row.names=FALSE)
colnames(feature.transform)
#the pls transformed 15 components 
pls_15components<-pls_feature$x.scores


#sparse PCA
#choose 10 or 20 features
library(elasticnet)
k=20
para_k=rep(0.5,k)
spca.feature<-spca(feature.transform,K=k,para=para_k,type="predictor",sparse="penalty")
#loadings The loadings of the sparse PCs
#pev Percentage of explained variance
#var.all Total variance of the predictors
dim(spca.feature$loadings)
spca.feature$pev
sum(spca.feature$pev)
library(nsprcomp)
spca.feature2<-nsprcomp(feature.transform, retx = TRUE, ncomp = 20)
dim(spca.feature2$rotation)
dim(spca.feature2$x)
prcomp.spca.20components<-spca.feature2$x
prcomp.spca.20features[1:5,]
colnames(prcomp.spca.20components)<-sapply(1:ncol(prcomp.spca.20components),function(x){gsub('PC','sPCA_Component',colnames(prcomp.spca.20components)[x])})
spca_20components_y<-cbind(prcomp.spca.20components,diseased_group)
write.csv(spca_20components_y, file = "sPCA_20components_diseasegroup.csv",row.names=FALSE)
##paired t test
#sparse PLS
library(spls)
#choose 10 or 20 features
#use cross validation to decide number of principal components
dat.resp.new=2-dat.resp
test.spls.cv<-cv.splsda(feature.transform,dat.resp.new, fold=10, K=10, eta=0.7, kappa=0.5,
                         classifier='LDA', scale.x=TRUE, plot.it=TRUE, n.core=6 )
dat.resp.new<-as.numeric(unlist(dat.resp.new))
f <- splsda( feature.transform,dat.resp.new, K=10, eta=0.5 )
pc.spls.10.features<-feature.transform%*%f$W
f.5 <- splsda( feature.transform,dat.resp.new, K=5, eta=0.5 )
pc.spls.5.features<-feature.transform%*%f.5$W

