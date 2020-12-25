library(cluster)
###Exercise 1
x1 = c("blue", T, T, T, 12)
x2 = c("red", F, F,NA, NA)
x3 = c("red", T, F,NA, 17)
x4 = c("green", T, F, F, 21)

####a
Gower.d.12=(1/1+1/1+1/1)/3
Gower.d.13=(1+0+1+0+(17-12)/9)/4
Gower.d.14=(1+0+1+1+1)/5
Gower.d.23=(0+1)/3
Gower.d.24=(1+1)/3
Gower.d.34=(1+0+(21-17)/9)/4

dist=c(Gower.d.12,Gower.d.13,Gower.d.23,Gower.d.14,Gower.d.24,Gower.d.34)
dist.mat <- matrix(0, nrow = 3, ncol = 3)
rownames(dist.mat)=c("x1","x2","x3")
colnames(dist.mat)=c("x2","x3","x4")
dist.mat[upper.tri(dist.mat, diag = TRUE)] <- dist
t(dist.mat)

###b
j.d.12=(1-0)
j.d.13=(1-1/2)
j.d.14=(1-1/2)
j.d.23=(1)
j.d.24=(1)
j.d.34=(1-1)

bGower.d.12=(1/1+3*j.d.12)/4
bGower.d.13=(1+3*j.d.13+(17-12)/9)/5
bGower.d.14=(1+3*j.d.14+1)/5
bGower.d.23=(0+3*j.d.23)/4
bGower.d.24=(1+3*j.d.24)/4
bGower.d.34=(1+3*j.d.34+(21-17)/9)/5

bdist=c(bGower.d.12,bGower.d.13,bGower.d.23, bGower.d.14, bGower.d.24, bGower.d.34)
bdist.mat <- matrix(0, nrow = 3, ncol = 3)
rownames(bdist.mat)=c("x1","x2","x3")
colnames(bdist.mat)=c("x2","x3","x4")
bdist.mat[upper.tri(bdist.mat, diag = TRUE)] <- bdist 
t(bdist.mat)



####c
d=as.data.frame(rbind(x1,x2,x3,x4))
d[,5]=as.numeric(as.character(d[,5]))
str(d)
daisy(d,metric = "gower")




#Exercise 3
library(factoextra)
library(cluster)
library(prabclus)
library(smacof)
data(tetragonula)

ta <- alleleconvert(strmatrix=tetragonula)
tai <- alleleinit(allelematrix=ta)

###########a
mdstai=mds(tai$distmat)


plot(mdstai$conf, xlab = "", ylab = "", asp = 1)

#Dendogram
library("ggplot2")
library("ggdendro")
library(ggpubr)

dist(mdstai$conf)
#####Complete
clusta <- hclust(dist(mdstai$conf),method="complete")
a=ggdendrogram(clusta, rotate = F, theme_dendro = F,labels = F)
print(a + ggtitle("Complete Linkage"))

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:30){
  tclusk[[k]] <- cutree(clusta,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=tai$distmat)
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=tai$distmat))$avg.width
}
plot(1:30,tasw,type="l",xlab="Number of clusters",ylab="ASW")
m=match(max(tasw[-1]),tasw[-1])+1
abline(v=m,col="red",lty="dashed")


clusta12=cutree(clusta,m)

clusta12f=as.factor(clusta12)

mdstai.df=as.data.frame(mdstai$conf)
mdstai.df= cbind(mdstai.df,clusta12f)

gg1=ggscatter( mdstai.df,
               x="D1", y="D2",
               color = "clusta12f", palette = "simpsons",
               shape = "clusta12f",
               ellipse = F, 
               mean.point = TRUE,
               star.plot = TRUE)
gg1+theme(legend.position = "none")




### Average
clustb <- hclust(dist(mdstai$conf),method="average")
b=ggdendrogram(clustb, rotate = F, theme_dendro = F,labels = F)
print(b + ggtitle("Average Linkage"))


tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:30){
  tclusk[[k]] <- cutree(clustb,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=dist(mdstai$conf))
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=dist(mdstai$conf)))$avg.width
}
plot(1:30,tasw,type="l",xlab="Number of clusters",ylab="ASW")
m=match(max(tasw[-1]),tasw[-1])+1
abline(v=m,col="red",lty="dashed")


clustb10=cutree(clustb,m)

clustb10f=as.factor(clustb10)

mdstai.df= cbind(mdstai.df,clustb10f)

gg1=ggscatter( mdstai.df,
               x="D1", y="D2",
               color = "clustb10f", palette = "simpsons",
               shape = "clustb10f",
               ellipse = F, 
               mean.point = TRUE,
               star.plot = TRUE)
gg1+theme(legend.position = "none")




#Ward's method
wclust <- hclust(dist(mdstai$conf),method="ward.D2")
w=ggdendrogram(wclust, rotate = F, theme_dendro = F,labels = F)
print(w + ggtitle("Ward’s method"))

tasw <- NA
tclusk <- list()
tsil <- list()
for (k in 2:30){
  tclusk[[k]] <- cutree(wclust,k)
  tsil[[k]] <- silhouette(tclusk[[k]],dist=dist(mdstai$conf))
  tasw[k] <- summary(silhouette(tclusk[[k]],dist=dist(mdstai$conf)))$avg.width
}
plot(1:30,tasw,type="l",xlab="Number of clusters",ylab="ASW")
m=match(max(tasw[-1]),tasw[-1])+1
abline(v=m,col="red",lty="dashed")


wclust9=cutree(wclust,m)

wclust9f=as.factor(wclust9)

mdstai.df= cbind(mdstai.df,wclust9f)

gg1=ggscatter( mdstai.df,
               x="D1", y="D2",
               color = "wclust9f", palette = "simpsons",
               shape = "wclust9f",
               ellipse = F, 
               mean.point = TRUE,
               star.plot = TRUE)
gg1+theme(legend.position = "none")


####PAM

pasw <- NA
pclusk <- list()
psil <- list()
for (k in 2:30){
  pclusk[[k]] <- pam(dist(mdstai$conf),k)
  # Computation of silhouettes:
  psil[[k]] <- silhouette(pclusk[[k]],dist=dist(mdstai$conf))
  # ASW needs to be extracted:
  pasw[k] <- summary(psil[[k]])$avg.width
}


plot(1:30,pasw,type="l",xlab="Number of clusters",ylab="ASW")
m=match(max(pasw[-1]),pasw[-1])+1
abline(v=m,col="red",lty="dashed")

pam7=as.factor(pclusk[[m]]$cluster)
mdstai.df= cbind(mdstai.df,pam7)
gg1=ggscatter( mdstai.df,
               x="D1", y="D2",
               color = "pam7", palette = "simpsons",
               shape = "pam7",
               ellipse = F, 
               mean.point = TRUE,
               star.plot = TRUE)
gg1+theme(legend.position = "none")

sil<- matrix(0, nrow = 1, ncol = 4)

colnames(sil)=c("Complete","Average","Ward","Pam")
rownames(sil)=c("ASW")
sil[1,]=c(0.4213881, 0.4213881,0.4245236,0.6401453)
sil
###########b
cg1 <- clusGap(mdstai$conf,kmeans,13,B=100,d.power=2,spaceH0="scaledPCA",nstart=100)
plot(cg1,main="")
a=print(cg1,method="globalSEmax",SE.factor=2)
b=print(cg1,method="Tibs2001SEmax",SE.factor=2)

kmeans=kmeans(mdstai$conf,7,nstart=100)
km=as.factor(kmeans$cluster)

mdstai.df= cbind(mdstai.df,km)
gg1=ggscatter( mdstai.df,
               x="D1", y="D2",
               color = "km", palette = "simpsons",
               shape = "km",
               ellipse = F, 
               mean.point = TRUE,
               star.plot = TRUE)
gg1+theme(legend.position = "none")+ggtitle("Kmeans 7 Clusters")
adjustedRandIndex(kmeans$cluster,pclusk[[7]]$cluster)

a=c(0.4213881, 0.4213881,0.4245236,0.6401453)