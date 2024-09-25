
library(vegan)
library(psych)
library(data.table)
rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}
setwd("D:/SC_bio")

gen=fread("SNP_TMP.txt",header=F)
POP=read.table("SC_group.txt",header=F)
LOCI=read.table("ID.txt",header=F)

#注释表头
rownames(gen)=as.character(POP$V1)
colnames(gen)=as.character(LOCI$V1)

env <- read.csv("SC_19climates.csv",head=T)
env$ID <- as.character(env$ID)
pred=env[,c(2:21)]

pk.rda <- rda(gen ~ ., data=pred, scale=T)
load.rda <- scores(pk.rda, choices=c(1:6), display="species")
head(load.rda)
write.table(load.rda,file="loaduniq_rda.txt", quote=F)
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)
  x[x < lims[1] | x > lims[2]]
}
cand1 <- outliers(load.rda[,1],3)
cand2 <- outliers(load.rda[,2],3)
cand3 <- outliers(load.rda[,3],3)
cand4 <- outliers(load.rda[,4],3)
cand5 <- outliers(load.rda[,5],3)
cand6 <- outliers(load.rda[,6],3)

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
cand4 <- cbind.data.frame(rep(4,times=length(cand4)), names(cand4), unname(cand4))
cand5 <- cbind.data.frame(rep(5,times=length(cand5)), names(cand5), unname(cand5))
cand6 <- cbind.data.frame(rep(6,times=length(cand6)), names(cand6), unname(cand6))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- colnames(cand4)<- colnames(cand5)<-colnames(cand6) <- c("axis","snp","loading")


write.table(cand1,file="result/cand1.txt",quote=F)
write.table(cand2,file="result/cand2.txt",quote=F)
write.table(cand3,file="result/cand3.txt",quote=F)
write.table(cand4,file="result/cand4.txt",quote=F)
write.table(cand5,file="result/cand5.txt",quote=F)
write.table(cand6,file="result/cand6.txt",quote=F)
ncand <- length(cand1) + length(cand2) + length(cand3)+length(cand4) + length(cand5) + length(cand6)

cand <- rbind(cand1, cand2, cand3,cand4, cand5, cand6)
cand$snp <- as.character(cand$snp)

# 初始化一个空的矩阵，用于存储相关性结果，并赋予列名

foo <- matrix(NA, nrow = nrow(cand), ncol = 20)

colnames(foo) <- c( "alt","bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")

gen <- as.data.frame(gen)

for (i in 1:nrow(cand)) {
  nam <- cand[i, "snp"]
  if (nam %in% colnames(gen)) {
    snp.gen <- gen[, nam]
    cor_values <- apply(pred, 2, function(x) cor(x, snp.gen))
    if (length(cor_values) == ncol(foo)) {
      foo[i, ] <- cor_values
    } else {
      warning(paste("相关性值的数量与 foo 矩阵的列数不匹配，位点：", nam))
    }
  } else {
    warning(paste("位点", nam, "在基因型数据中不存在，将跳过此位点。"))
  }
}

cand <- cbind.data.frame(cand,foo)
cand <- cand[!duplicated(cand$snp),]

result <- data.frame(matrix(ncol = 2, nrow = nrow(cand)))
colnames(result) <- c("predictor", "correlation")

# 循环计算每个基因位点与环境因子的相关性
for (i in 1:nrow(cand)) {
  bar <- cand[i, 4:ncol(cand)]
  result[i, "predictor"] <- names(which.max(abs(bar)))   # 获取最相关的环境因子列名
  result[i, "correlation"] <- max(abs(bar))             # 获取相关性的最大值
}

write.table(cand,file="result/cand.txt",quote=F)
write.table(result, file = "result/cand_result.txt", quote = FALSE, row.names = FALSE)

pdf("RDA.pdf")
opar<-par(no.readonly=TRUE)
par(mfrow=c(3,2))
hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3")
hist(load.rda[,4], main="Loadings on RDA4")
hist(load.rda[,5], main="Loadings on RDA5")
hist(load.rda[,6], main="Loadings on RDA6")
par(opar)
dev.off()
