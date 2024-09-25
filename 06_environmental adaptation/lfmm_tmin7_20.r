#!/usr/bin/env Rscript
library(LEA)
library(qvalue)
#library(tidyverse)
setwd("/work/OrthoTree/lfmm")
obj.lfmm = lfmm("EN_chr20.lfmm", "./env/tmin7.env", K = 3, rep = 5,CPU = 64, project="new")
zs = z.scores(obj.lfmm, K = 3)
zs.median = apply(zs, MARGIN = 1, median)
lambda = median(zs.median^2)/qchisq(0.5, df = 1)
lambda
adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE)

tiff("tmin7.chr20.pvalue.tiff")
hist = hist(adj.p.values, col = "red")
plot(hist)
dev.off() 

L = 223309
q = 0.1
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates.bh = order(adj.p.values)[w]

tiff("tmin7.chr20.qvalue.tiff")
plot(qvalue(adj.p.values))
dev.off() 

candidates.qv = which(qvalue(adj.p.values, fdr = .1)$signif)

write.table(candidates.bh,file="./tmin7.chr20.pvalue.candiloci.txt")
write.table(candidates.qv,file="./tmin7.chr20.qvalue.candiloci.txt")
write.table(adj.p.values,file="./tmin7.chr20.pvalue.txt")
